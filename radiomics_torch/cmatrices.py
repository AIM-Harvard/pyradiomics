from typing import List, Optional
import torch
from itertools import product
import torch.nn.functional as F

### GLCM

def get_angle_count_torch(size, distances, Nd, bidirectional=False, force2Ddim=-1):
    """
    逐行翻译自 C 版 get_angle_count
    size: 序列[int]，长度 Nd
    distances: 序列[int]
    Nd: 维度数
    bidirectional: False 对应 C 里的 0（只要单向），True 对应 1（双向）
    force2Ddim: 被排除的 out-of-plane 维度，或 -1（不排除）
    """
    size = [int(s) for s in size]
    distances = [int(d) for d in distances]
    Na = 0

    for d in distances:
        if d < 1:
            return 0

        Na_d = 1
        Na_dd = 1
        for dim_idx in range(Nd):
            if dim_idx == force2Ddim:
                continue

            if d < size[dim_idx]:
                Na_d *= (2 * d + 1)
                Na_dd *= (2 * d - 1)
            else:
                Na_d  *= (2 * (size[dim_idx] - 1) + 1)
                Na_dd *= (2 * (size[dim_idx] - 1) + 1)

        Na += (Na_d - Na_dd)

    if not bidirectional:
        Na //= 2

    return Na


def build_angles_torch(size, distances, Nd, Ndist, force2Ddim, Na):
    """
    逐行翻译自 C 版 build_angles
    返回 angles: list[list[int]]，形状 (Na, Nd)
    """
    size = [int(s) for s in size]
    distances = [int(d) for d in distances]

    # 找到 max_distance
    max_distance = 0
    for d in distances:
        if d < 1:
            raise ValueError("Invalid distance encountered (<1)")
        if max_distance < d:
            max_distance = d

    n_offsets = 2 * max_distance + 1

    # offset_stride，和 C 里一样
    offset_stride = [0] * Nd
    offset_stride[Nd - 1] = 1
    for dim_idx in range(Nd - 2, -1, -1):
        offset_stride[dim_idx] = offset_stride[dim_idx + 1] * n_offsets

    angles = [[0] * Nd for _ in range(Na)]

    new_a_idx = 0  # 控制 offset 组合
    a_idx = 0      # 真正写入到 angles 的 index（只对有效角度 ++）

    while a_idx < Na:
        a_dist = 0

        for dim_idx in range(Nd):
            offset = max_distance - (new_a_idx // offset_stride[dim_idx]) % n_offsets

            if ((dim_idx == force2Ddim and offset != 0) or
                offset >= size[dim_idx] or
                offset <= -size[dim_idx]):
                a_dist = -1    # invalid angle
                break

            angles[a_idx][dim_idx] = offset

            if a_dist < offset:
                a_dist = offset
            elif a_dist < -offset:
                a_dist = -offset

        new_a_idx += 1

        # a_dist < 1：无效（<0 不合法 or ==0 是 (0,0,0)）
        if a_dist < 1:
            continue

        # 检查这个距离是不是在 distances 里面
        ok = False
        for d in distances:
            if a_dist == d:
                ok = True
                break

        if ok:
            a_idx += 1

    return angles


def build_angles_arr_torch(distances_obj, size, Nd, force2Ddimension, bidirectional=False):
    """
    高层封装，对应 C 版 build_angles_arr

    distances_obj: None 或 1D 序列
    size: 图像 size，len=size Nd
    Nd: 维数
    force2Ddimension: 被排除的维度 index 或 -1
    bidirectional: False 单向（默认，和 PyRadiomics 一致），True 双向

    返回:
      angles: LongTensor (Na, Nd)
      Na: int
    """
    if distances_obj is None:
        distances = [1]
    else:
        distances = [int(d) for d in distances_obj]

    Ndist = len(distances)
    if Ndist == 0:
        raise ValueError("distances must be non-empty")

    Na = get_angle_count_torch(size, distances, Nd, bidirectional, force2Ddimension)
    if Na == 0:
        raise RuntimeError("Error getting angle count. Check distances and image size.")

    angles_list = build_angles_torch(size, distances, Nd, Ndist, force2Ddimension, Na)
    angles = torch.tensor(angles_list, dtype=torch.long)
    return angles, Na

def set_bb_torch(v, size, voxels, Nd, Nvox, kernelRadius, force2Ddimension):
    """
    Python/Torch version of set_bb (cmatrices.c)

    size: 序列[int]，len Nd
    voxels: None 或 LongTensor (Nd, Nvox)
    返回:
      bb_lo, bb_hi: 各是长度 Nd 的 list[int]，闭区间 [lo, hi]
    """
    size = [int(s) for s in size]
    bb_lo = [0] * Nd
    bb_hi = [s - 1 for s in size]

    if voxels is not None:
        # voxels 形状 (Nd, Nvox)
        for d in range(Nd):
            if d == force2Ddimension:
                coord = int(voxels[d, v])
                bb_lo[d] = coord
                bb_hi[d] = coord
            else:
                coord = int(voxels[d, v])
                lo = coord - kernelRadius
                hi = coord + kernelRadius
                if lo < 0:
                    lo = 0
                if hi >= size[d]:
                    hi = size[d] - 1
                bb_lo[d] = lo
                bb_hi[d] = hi

    return bb_lo, bb_hi

def _build_slices_for_offset_torch(shape, offset):
    """
    给定 shape=(N0,...,Nd-1) 和 offset=(o0,...,od-1)，
    生成 src / dst slice 列表，对应 C 里 cur_idx + angles 检查边界的效果。
    """
    src_slices = []
    dst_slices = []
    for n, off in zip(shape, offset):
        off = int(off)
        if off > 0:
            src_slices.append(slice(0, n - off))
            dst_slices.append(slice(off, n))
        elif off < 0:
            src_slices.append(slice(-off, n))
            dst_slices.append(slice(0, n + off))
        else:
            src_slices.append(slice(None))
            dst_slices.append(slice(None))
    return tuple(src_slices), tuple(dst_slices)


def _calculate_glcm_local_torch(image, mask, angles, Ng):
    """
    在一个给定 bounding box 内计算 GLCM，
    是 C 版 `calculate_glcm` 的矢量化等价实现。

    image, mask: 相同 shape 的 tensor
    angles: (Na, Nd) LongTensor
    Ng: 灰度级数 (image 值应该在 1..Ng, 0/<=0 会被忽略)

    返回 glcm: (Ng, Ng, Na)
    """
    if not torch.is_tensor(image):
        image = torch.as_tensor(image)
    if not torch.is_tensor(mask):
        mask = torch.as_tensor(mask)

    device = image.device
    image = image.to(device=device, dtype=torch.int64)
    mask = mask.to(device=device, dtype=torch.bool)
    angles = angles.to(device=device, dtype=torch.long)

    assert image.shape == mask.shape
    Nd = image.dim()
    Na = angles.shape[0]
    assert angles.shape[1] == Nd

    shape = image.shape
    glcm = torch.zeros((Ng, Ng, Na), device=device, dtype=torch.float64)

    for a in range(Na):
        offset = angles[a]
        src_slices, dst_slices = _build_slices_for_offset_torch(shape, offset)

        src_mask = mask[src_slices]
        dst_mask = mask[dst_slices]
        valid = src_mask & dst_mask
        if not torch.any(valid):
            continue

        src_vals = image[src_slices][valid]
        dst_vals = image[dst_slices][valid]

        good = (
            (src_vals > 0)
            & (dst_vals > 0)
            & (src_vals <= Ng)
            & (dst_vals <= Ng)
        )
        if not torch.any(good):
            continue

        i_vals = src_vals[good] - 1
        j_vals = dst_vals[good] - 1

        idx = i_vals * Ng + j_vals  # flatten 到 [0, Ng*Ng)
        bc = torch.bincount(idx, minlength=Ng * Ng).to(dtype=torch.float64)
        glcm[..., a] = bc.view(Ng, Ng)

    return glcm

def calculate_glcm_torch(
    image,
    mask,
    distances,
    Ng,
    force2D,
    force2Ddimension,
    kernelRadius=0,
    voxels=None,
):
    """
    高层 PyTorch 版，相当于 C 扩展的 cmatrices_calculate_glcm + calculate_glcm.

    参数说明基本和 C 版一致：
      image : Tensor / array-like，已经离散化（1..Ng），0/<=0 当作无效
      mask  : Tensor / array-like，非零为 ROI
      distances : None 或 [d1, d2, ...]
      Ng    : int
      force2D : bool/int，False -> 不强制 2D，True -> 用 force2Ddimension
      force2Ddimension : int，被排除的 out-of-plane 维度（如果 force2D=True）
      kernelRadius : int, voxel-based 半径
      voxels : None 或 shape (Nd, Nvox) 的 voxel indices（列向量）

    返回:
      glcm_all : (Nvox, Ng, Ng, Na) float64, 在 image.device 上
      angles   : (Na, Nd) int64, 在 CPU 上（行为上接近 numpy）
    """
    # 转 tensor，推断 device
    if not torch.is_tensor(image):
        image = torch.as_tensor(image)
    if not torch.is_tensor(mask):
        mask = torch.as_tensor(mask)

    if image.shape != mask.shape:
        raise ValueError("Image and mask must have the same shape")

    device = image.device
    Nd = image.dim()
    size = list(image.shape)

    # force2D 逻辑和 C 版一样：如果不 force2D，就把维度设为 -1
    if not force2D:
        f2d_dim = -1
    else:
        f2d_dim = int(force2Ddimension)

    # distances 处理
    if distances is None:
        distances_list = [1]
    else:
        distances_list = [int(d) for d in distances]
    if len(distances_list) == 0:
        raise ValueError("distances must be non-empty")

    # 生成 angles（等价 build_angles_arr）
    angles_cpu, Na = build_angles_arr_torch(
        distances_list, size, Nd, f2d_dim, bidirectional=False
    )
    # 用在计算时 copy 一份到 device；返回给上层时保留 CPU 版
    angles = angles_cpu.to(device=device)

    # 处理 voxel-based / segment-based
    voxels_t = None
    if voxels is not None:
        voxels_t = torch.as_tensor(voxels, dtype=torch.long)
        if voxels_t.dim() != 2 or voxels_t.shape[0] != Nd:
            raise ValueError("voxels must have shape (Nd, Nvox)")
        Nvox = int(voxels_t.shape[1])
    else:
        Nvox = 1

    # 分配输出：和 C 版一样 (Nvox, Ng, Ng, Na)
    glcm_all = torch.zeros((Nvox, Ng, Ng, Na), device=device, dtype=torch.float64)

    # segment-based: Nvox=1 且 bb=全图
    # voxel-based: 按 set_bb 的逻辑为每个 v 建一个局部 bb，然后在 bb 上跑 _calculate_glcm_local_torch
    for v in range(Nvox):
        bb_lo, bb_hi = set_bb_torch(
            v, size, voxels_t, Nd, Nvox, kernelRadius, f2d_dim
        )
        slices = tuple(slice(bb_lo[d], bb_hi[d] + 1) for d in range(Nd))
        sub_image = image[slices]
        sub_mask = mask[slices]

        glcm_local = _calculate_glcm_local_torch(sub_image, sub_mask, angles, Ng)
        glcm_all[v] = glcm_local

    return glcm_all, angles_cpu

### gldm

def calculate_gldm_torch(
    image: torch.Tensor,
    mask: torch.Tensor,
    distances,
    Ng: int,
    alpha: int,
    force2D: bool,
    force2Ddimension: int,
    kernelRadius: int = 0,
    voxels: torch.Tensor | None = None,
    device=None,
    dtype=torch.float32,
):
    """
    Torch 版 GLDM 计算，逻辑对应 cmatrices_calculate_gldm + calculate_gldm.

    参数：
      image: int Tensor, 形状 (N1, N2, ..., Nd)，值范围 [1..Ng]
      mask:  bool/0-1 Tensor, 同 shape
      distances: 与 C 版 distances_obj 对应（list/1D tensor）
      Ng:   灰度级数
      alpha: 依赖阈值 (|i-j|<=alpha 视为 dependent)
      force2D, force2Ddimension: 与 C 版含义一致
      kernelRadius: voxel-based 时局部核半径，>0 才启用 voxels
      voxels: None 或 LongTensor (Nd, Nvox)，与 C 版 try_parse_voxels_arr 一致

    返回：
      gldm: Tensor, shape (Nvox, Ng, Na*2+1), dtype=dtype
    """
    if device is None:
        device = image.device

    image = image.to(device)
    mask = mask.to(device)

    # Nd, size
    Nd = image.ndim
    size = list(image.shape)

    # 如果不强制 2D，则设为 -1（和 C 版一样）
    if not force2D:
        force2Ddimension = -1

    # ---- 构建 angles（对应 build_angles_arr, 最后一个参数 C 里是 1：双向）----
    angles, Na = build_angles_arr_torch(
        distances_obj=distances,
        size=size,
        Nd=Nd,
        force2Ddimension=force2Ddimension,
        bidirectional=True,   # C 里最后一个参数是 1
    )  # angles: (Na, Nd), long

    # ---- Nvox / voxels 形状 ----
    if voxels is not None:
        voxels = voxels.to(device)
        # C 版 expect 形状 Nd x Nvox
        assert voxels.shape[0] == Nd, f"voxels must be (Nd, Nvox), got {voxels.shape}"
        Nvox = voxels.shape[1]
    else:
        Nvox = 1  # ROI-based, 整块 ROI 一次一个 GLDM

    # ---- gldm 输出 shape ----
    max_dep = Na * 2 + 1
    gldm = torch.zeros((Nvox, Ng, max_dep), dtype=dtype, device=device)

    # ---- 构造 strides (元素为单位，不是字节) ----
    # C 版 try_parse_arrays 里算的是 "元素 stride"
    strides = [0] * Nd
    run = 1
    for d in range(Nd - 1, -1, -1):
        strides[d] = run
        run *= size[d]

    # ---- flatten image / mask 方便用线性 index ----
    image_flat = image.reshape(-1).to(torch.int64)
    mask_flat = mask.reshape(-1).to(torch.bool)

    # ---- 主循环：每个 voxel / ROI 一个 GLDM ----
    for v in range(Nvox):
        # bb 的定义与 C 版 set_bb 完全一致（你已经在 set_bb_torch 里实现了）
        bb_lo, bb_hi = set_bb_torch(
            v=v,
            size=size,
            voxels=voxels,              # None 表示用整个 ROI
            Nd=Nd,
            Nvox=Nvox,
            kernelRadius=kernelRadius,
            force2Ddimension=force2Ddimension,
        )
        gldm_v = _calculate_gldm_single_torch(
            image_flat=image_flat,
            mask_flat=mask_flat,
            size=size,
            bb_lo=bb_lo,
            bb_hi=bb_hi,
            strides=strides,
            angles=angles,
            Na=Na,
            Nd=Nd,
            Ng=Ng,
            alpha=alpha,
            max_dep=max_dep,
            dtype=dtype,
            device=device,
        )
        gldm[v] = gldm_v

    return gldm

def _calculate_gldm_single_torch(
    image_flat: torch.Tensor,
    mask_flat: torch.Tensor,
    size,
    bb_lo,
    bb_hi,
    strides,
    angles: torch.Tensor,  # (Na, Nd)
    Na: int,
    Nd: int,
    Ng: int,
    alpha: int,
    max_dep: int,
    dtype=torch.float32,
    device=None,
):
    """
    对应 C 版 calculate_gldm，但只负责当前一个 bb / voxel 的 gldm 计算。

    返回:
      gldm_v: (Ng, max_dep)
    """
    if device is None:
        device = image_flat.device

    gldm_v = torch.zeros((Ng, max_dep), dtype=dtype, device=device)

    # 构造每个维度上的 index 范围
    ranges = [range(bb_lo[d], bb_hi[d] + 1) for d in range(Nd)]

    # 遍历 bounding box 内所有 voxel（对应 C 版 i 从 bb 下界扫到上界）
    for idx in product(*ranges):
        # 线性 index i = sum(idx[d] * strides[d])
        i_lin = 0
        for d in range(Nd):
            i_lin += idx[d] * strides[d]

        if not mask_flat[i_lin]:
            continue

        center_val = int(image_flat[i_lin].item())
        if center_val <= 0:
            # 对应 C 版: if (image[i] <= 0 || gldm_idx >= gldm_idx_max) return 0;
            # 这里我们严格一些，直接跳过或抛异常都可以：
            # raise RuntimeError("GLDM: image value <= 0 encountered.")
            continue

        dep = 0  # 依赖邻居数

        # 遍历所有角度 a（对应 a=0..Na-1）
        for a in range(Na):
            # 当前角 offset 各维
            neighbor_idx = [0] * Nd
            out_of_range = False
            for d in range(Nd):
                off = int(angles[a, d].item())
                coord = idx[d] + off
                # C 版：如果超出 bb，则 j 置为 i 并 break
                if coord < bb_lo[d] or coord > bb_hi[d]:
                    out_of_range = True
                    break
                neighbor_idx[d] = coord

            if out_of_range:
                continue

            # 线性 index j
            j_lin = 0
            for d in range(Nd):
                j_lin += neighbor_idx[d] * strides[d]

            # j==i 的情况在 C 里表示“没有移动”，这里不会出现（因为只要有一个 offset 非零）
            if not mask_flat[j_lin]:
                continue

            diff = abs(int(image_flat[i_lin].item()) - int(image_flat[j_lin].item()))
            if diff <= alpha:
                dep += 1

        # 写入 gldm，C 版：
        #   gldm_idx = dep + (image[i]-1) * (Na * 2 + 1);
        #   gldm[gldm_idx]++;
        if dep < 0 or dep >= max_dep:
            # 理论上 dep in [0, Na*2]，这里加个保护
            continue

        gray_idx = center_val - 1  # 0-based
        if gray_idx < 0 or gray_idx >= Ng:
            # out of range, 对应 C 版直接 return 0
            # 可以选择 raise/continue，这里为了鲁棒先跳过
            continue

        gldm_v[gray_idx, dep] += 1.0

    return gldm_v


### glrlm

import torch
import torch.nn.functional as F


def shift3d(x, direction):
    """
    x: (B=1, C=1, Z, Y, X)
    direction: (dz, dy, dx)
    """
    dz, dy, dx = direction
    _, _, Z, Y, X = x.shape

    pad = [
        max(dx,0), max(-dx,0),   # left/right
        max(dy,0), max(-dy,0),   # top/bottom
        max(dz,0), max(-dz,0),   # front/back
    ]

    x_pad = F.pad(x, pad)
    z0 = max(-dz,0)
    y0 = max(-dy,0)
    x0 = max(-dx,0)

    return x_pad[:, :, z0:z0+Z, y0:y0+Y, x0:x0+X]


def glrlm_vectorized_3d(
    image,     # (Z,Y,X), int64, values ∈ [1..Ng]
    mask,      # (Z,Y,X), bool
    Ng: int,
    Nr: int,
    directions,   # list of (dz,dy,dx)
    device=None,
):
    """
    返回 P_glrlm: (Ng, Nr, Na) —— 不含 batch 维。
    """
    if device is None:
        device = image.device

    image = image.to(device=device, dtype=torch.int64)
    mask  = mask.to(device=device, dtype=torch.bool)

    Z, Y, X = image.shape
    Na = len(directions)

    # output
    P = torch.zeros((Ng, Nr, Na), dtype=torch.float32, device=device)

    # (1,1,Z,Y,X)
    img5 = image.unsqueeze(0).unsqueeze(0)
    msk5 = mask.unsqueeze(0).unsqueeze(0)

    # flatten helpers
    def flatten_runs(run_continue, vals):
        """
        run_continue: (N,) bool
        vals:         (N,) int64
        返回 run_start_idx, run_len, gray
        """
        N = run_continue.numel()
        rc = run_continue
        rc_prev = torch.zeros_like(rc)
        rc_prev[1:] = rc[:-1]

        run_start = rc & (~rc_prev)
        run_end   = rc_prev & (~rc)

        idx = torch.arange(N, device=device)

        start_idx = idx[run_start]
        end_idx   = idx[run_end]

        # 如果最后 run 没有结束，补一个结束点
        if end_idx.numel() < start_idx.numel():
            end_idx = torch.cat([end_idx, idx[-1:].clone()])

        run_lengths = (end_idx - start_idx + 1).clamp(max=Nr).long()
        run_gray = (vals[start_idx] - 1).clamp(0, Ng-1).long()

        return run_gray, run_lengths

    # ---- 对每个方向向量化计算 ----
    for a, (dz,dy,dx) in enumerate(directions):
        img_s = shift3d(img5, (dz,dy,dx))[0,0]  # (Z,Y,X)
        msk_s = shift3d(msk5, (dz,dy,dx))[0,0]

        same = (image == img_s)
        valid = mask & msk_s
        run_continue = (same & valid).flatten()
        vals = image.flatten()

        if run_continue.any():
            run_gray, run_lengths = flatten_runs(run_continue, vals)

            # bincount accumulate
            idx = run_gray * Nr + (run_lengths - 1)
            bc = torch.bincount(idx, minlength=Ng*Nr)

            P[:,:,a] = bc.view(Ng, Nr)

    return P


def calculate_glrlm_torch(
    image: torch.Tensor,
    mask: torch.Tensor,
    Ng: int,
    Nr: int,
    force2D: bool,
    force2Ddimension: int,
    kernelRadius: int = 0,
    voxels: torch.Tensor | None = None,
    dtype=torch.float32,
):
    """
    Torch 版 GLRLM 计算，对应 C 函数 cmatrices_calculate_glrlm + calculate_glrlm.

    参数
    ----
    image : int Tensor, 形状 (N1, N2, ..., Nd)，值范围 [1..Ng]
    mask  : uint8/bool Tensor, same shape as image
    Ng    : 灰度级数
    Nr    : 运行长度最大值（run length bins 个数）
    force2D : 是否强制 2D
    force2Ddimension : 被忽略的维度 index 或 -1
    kernelRadius : voxel-based 提取时的核半径（>0 时使用 voxels）
    voxels : None 或 LongTensor (Nd, Nvox)，和 C 版 try_parse_voxels_arr 约定一致
    device : 计算设备
    dtype  : glrlm 的 dtype，一般用 torch.double

    返回
    ----
    glrlm : Tensor, 形状 (Nvox, Ng, Nr, Na)
    angles: LongTensor, 形状 (Na, Nd)
    """
    device = image.device

    Nd = image.ndim
    size = list(image.shape)

    # 对齐 C 逻辑：如果不强制 2D，就把 force2Ddimension 设为 -1
    if not force2D:
        force2Ddimension = -1

    # ---- 生成角度数组，对应 C 里的 build_angles_arr(..., bidirectional=0) ----
    angles, Na = build_angles_arr_torch(
        distances_obj=None,       # GLRLM 固定距离
        size=size,
        Nd=Nd,
        force2Ddimension=force2Ddimension,
        bidirectional=False,      # C 里最后一个参数是 0
    )  # angles: (Na, Nd), long

    # ---- voxel-based / ROI-based ----
    if voxels is not None:
        voxels = voxels.to(device=device, dtype=torch.int64)
        assert voxels.shape[0] == Nd, f"voxels must be (Nd, Nvox), got {voxels.shape}"
        Nvox = voxels.shape[1]
    else:
        Nvox = 1
        voxels = None

    # ---- 构造 element-strides，和 C 版 try_parse_arrays 保持一致 ----
    strides = [0] * Nd
    run = 1
    for d in range(Nd - 1, -1, -1):
        strides[d] = run
        run *= size[d]
    Ni = run  # 整个数组的元素个数

    # ---- flatten image / mask ----
    image_flat = image.reshape(-1)
    mask_flat = mask.reshape(-1)

    # ---- 输出 glrlm ----
    glrlm = torch.zeros((Nvox, Ng, Nr, Na), dtype=dtype, device=device)

    bb_cache: dict[tuple, torch.Tensor] = {}

    for v in range(Nvox):
        # 1) 计算当前 voxel 的 bounding box
        bb_lo, bb_hi = set_bb_torch(
            v=v,
            size=size,
            voxels=voxels,
            Nd=Nd,
            Nvox=Nvox,
            kernelRadius=kernelRadius,
            force2Ddimension=force2Ddimension,
        )  # bb_lo[d], bb_hi[d] 是闭区间

        # 2) 用 (bb_lo, bb_hi) 作为 key 做缓存
        key = tuple(bb_lo + bb_hi)  # 例如 (z0,y0,x0,z1,y1,x1)

        if key in bb_cache:
            glrlm_v = bb_cache[key]
        else:
            glrlm_v = _calculate_glrlm_single_torch(
                image_flat=image_flat,
                mask_flat=mask_flat,
                size=size,
                bb_lo=bb_lo,
                bb_hi=bb_hi,
                strides=strides,
                angles=angles,
                Na=Na,
                Nd=Nd,
                Ng=Ng,
                Nr=Nr,
                Ni=Ni,
                dtype=dtype,
            )
            bb_cache[key] = glrlm_v

        glrlm[v] = glrlm_v

    return glrlm, angles


def _calculate_glrlm_single_torch(
    image_flat: torch.Tensor,
    mask_flat: torch.Tensor,
    size,
    bb_lo,
    bb_hi,
    strides,
    angles: torch.Tensor,  # (Na, Nd)
    Na: int,
    Nd: int,
    Ng: int,
    Nr: int,
    Ni: int,
    dtype=torch.float32,
):
    """
    对应 C 版 calculate_glrlm 的单 voxel/ROI 版本。
    返回 glrlm_v: (Ng, Nr, Na)
    """
    device = image_flat.device

    glrlm_v = torch.zeros((Ng, Nr, Na), dtype=dtype, device=device)

    # 计算 start_i：bounding box 左下角在展平数组上的 index
    start_i = 0
    for d in range(Nd):
        start_i += bb_lo[d] * strides[d]

    glrlm_idx_max = Ng * Nr * Na

    # 为了方便访问角度偏移，先搬到 CPU 的普通 Python int
    angles_np = angles.to("cpu").numpy()

    # ---- 遍历每个 angle ----
    for a in range(Na):
        # 找出当前角度的“moving dimensions”和起始 index
        mDims = []       # 移动维度的 index
        mDim_start = []  # 对应维度的起始坐标（bb_lo 或 bb_hi）
        for d in range(Nd):
            off = int(angles_np[a, d])
            if off != 0:
                if off > 0:
                    mDim_start.append(bb_lo[d])
                else:
                    mDim_start.append(bb_hi[d])
                mDims.append(d)
        cnt_mDim = len(mDims)

        multiElement = False  # 是否存在 run length > 1 的 run（对应 C 的 multiElement）

        # ---- 遍历 bounding box 范围内的所有线起点候选 i ----
        i = start_i
        while i < Ni:
            # 先把 i 调整到 bounding box 内（除维度 0 以外）
            for d in range(Nd - 1, 1, -1):  # 从 Nd-1 到 2
                cur_idx = (i % strides[d - 1]) // strides[d]
                if cur_idx > bb_hi[d]:
                    i += (size[d] - cur_idx + bb_lo[d]) * strides[d]
                elif cur_idx < bb_lo[d]:
                    i += (bb_lo[d] - cur_idx) * strides[d]

            if Nd > 1:
                # 单独处理 d = 1
                cur_idx = (i % strides[0]) // strides[1]
                if cur_idx > bb_hi[1]:
                    i += (size[1] - cur_idx + bb_lo[1]) * strides[1]
                elif cur_idx < bb_lo[1]:
                    i += (bb_lo[1] - cur_idx) * strides[1]

            cur_idx0 = i // strides[0]
            if cur_idx0 > bb_hi[0]:
                break  # 超出 bb，结束

            # 判断 i 是否是任何 moving dimension 上的“起始 voxel”
            start_voxel_valid = False
            for md, d in enumerate(mDims):
                if d == 0:
                    cur_idx = i // strides[d]
                else:
                    cur_idx = (i % strides[d - 1]) // strides[d]
                if cur_idx == mDim_start[md]:
                    start_voxel_valid = True
                    break

            if not start_voxel_valid:
                # 按 C 的逻辑，把 i 跳到最后一个 moving dimension 的下一个有效 start 位置
                if cnt_mDim == 0:
                    # 这个角度没有移动维度（理论上不该发生），跳出
                    break
                md = cnt_mDim - 1
                d = mDims[md]

                if d == 0:
                    cur_idx = i // strides[d]
                else:
                    cur_idx = (i % strides[d - 1]) // strides[d]

                delta = (mDim_start[md] - cur_idx + size[d]) % size[d]
                i += delta * strides[d]

                # 更新低维度，确保不出 bb
                if d > 1:
                    d -= 1
                    while d > 0:
                        cur_idx = (i % strides[d - 1]) // strides[d]
                        if cur_idx > bb_hi[d]:
                            i += (size[d] - cur_idx + bb_lo[d]) * strides[d]
                        elif cur_idx < bb_lo[d]:
                            i += (bb_lo[d] - cur_idx) * strides[d]
                        d -= 1

                if i // strides[0] > bb_hi[0]:
                    break
                # 跳完后继续 while i 循环
                continue

            # ---- 真正开始沿 angle 方向跑 run ----
            j = i
            gl = -1
            rl = 0
            elements = 0

            # do ... while (j != i) 结构
            first_iter = True
            while first_iter or (j != i):
                first_iter = False

                if mask_flat[j]:
                    elements += 1
                    cur_val = int(image_flat[j].item())
                    if gl == -1:
                        gl = cur_val
                        rl = 0
                    elif cur_val == gl:
                        rl += 1
                    else:
                        # 结束上一段 run，写入 glrlm
                        glrlm_idx = a + rl * Na + (gl - 1) * Na * Nr
                        if gl <= 0 or glrlm_idx >= glrlm_idx_max:
                            raise RuntimeError("GLRLM index out of range (run break)")
                        # 把一维 index转为 (gl_idx, rl_idx, a_idx)
                        gl_idx = gl - 1
                        rl_idx = rl
                        glrlm_v[gl_idx, rl_idx, a] += 1.0

                        gl = cur_val
                        rl = 0
                elif gl > -1:
                    # 当前 mask 为 0，结束 run
                    glrlm_idx = a + rl * Na + (gl - 1) * Na * Nr
                    if gl <= 0 or glrlm_idx >= glrlm_idx_max:
                        raise RuntimeError("GLRLM index out of range (mask break)")
                    gl_idx = gl - 1
                    rl_idx = rl
                    glrlm_v[gl_idx, rl_idx, a] += 1.0
                    gl = -1
                    rl = 0

                # 前进一步
                for md, d in enumerate(mDims):
                    if d == 0:
                        cur_idx = j // strides[d]
                    else:
                        cur_idx = (j % strides[d - 1]) // strides[d]

                    step = int(angles_np[a, d])
                    if cur_idx + step < bb_lo[d] or cur_idx + step > bb_hi[d]:
                        # 出边界，标志结束
                        j = i
                        break
                    j += step * strides[d]
                # while 结束条件由 (j != i) 控制

            # run 结束后如果 gl 还有效，需要再写一次
            if gl > -1:
                glrlm_idx = a + rl * Na + (gl - 1) * Na * Nr
                if gl <= 0 or glrlm_idx >= glrlm_idx_max:
                    raise RuntimeError("GLRLM index out of range (final run)")
                gl_idx = gl - 1
                rl_idx = rl
                glrlm_v[gl_idx, rl_idx, a] += 1.0

            if elements > 1:
                multiElement = True

            # 下一个 i
            i += 1

        # 处理“该 angle 只有 run length 1”的情况
        if not multiElement:
            # 对应 C 代码：
            # for (glrlm_idx = 0; glrlm_idx < Ng; glrlm_idx++)
            #   glrlm[glrlm_idx * Nr * Na + a] = 0;
            # 即，将 runlength index 0 的所有 gray level 在该 angle 上置 0
            glrlm_v[:, 0, a] = 0.0

    return glrlm_v



### First Order

def _get_angle_count_core(size_list, dist_list, bidirectional=True, force2Ddim=-1):
    """和 get_angle_count_np 完全同逻辑，只用 Python 标量和 list。"""
    Nd = len(size_list)
    Ndist = len(dist_list)

    Na = 0
    for dist_idx in range(Ndist):
        d = int(dist_list[dist_idx])
        if d < 1:
            return 0

        Na_d = 1
        Na_dd = 1
        for dim_idx in range(Nd):
            if dim_idx == force2Ddim:
                continue

            if d < size_list[dim_idx]:
                Na_d *= (2 * d + 1)
                Na_dd *= (2 * d - 1)
            else:
                max_step = 2 * (size_list[dim_idx] - 1) + 1
                Na_d *= max_step
                Na_dd *= max_step

        Na += (Na_d - Na_dd)

    if not bidirectional:
        Na //= 2

    return Na


def _build_angles_core_to_tensor(size_list, dist_list, Nd, Ndist, force2Ddim, Na, device, dtype):
    """和 build_angles_np 完全同逻辑，但写入的是 Torch Tensor。"""
    # 找最大距离
    max_distance = 0
    for d in dist_list:
        if d < 1:
            raise ValueError("Invalid distance (< 1) encountered in build_angles")
        if max_distance < d:
            max_distance = d

    n_offsets = 2 * max_distance + 1

    # offset_stride[Nd-1] = 1; 其余反向累乘
    offset_stride = [0] * Nd
    offset_stride[Nd - 1] = 1
    for dim_idx in range(Nd - 2, -1, -1):
        offset_stride[dim_idx] = offset_stride[dim_idx + 1] * n_offsets

    angles = torch.empty((Na, Nd), dtype=dtype, device=device)

    new_a_idx = 0
    a_idx = 0

    while a_idx < Na:
        a_dist = 0
        valid_angle = True

        for dim_idx in range(Nd):
            offset = max_distance - ((new_a_idx // offset_stride[dim_idx]) % n_offsets)

            if ((dim_idx == force2Ddim and offset != 0) or
                offset >= size_list[dim_idx] or
                offset <= -size_list[dim_idx]):
                a_dist = -1
                valid_angle = False
                break

            angles[a_idx, dim_idx] = int(offset)

            if a_dist < offset:
                a_dist = offset
            elif a_dist < -offset:
                a_dist = -offset

        new_a_idx += 1

        if a_dist < 1:
            continue
        if not valid_angle:
            continue

        for d in dist_list:
            if a_dist == d:
                a_idx += 1
                break

    return angles


def generate_angles_torch(
    size,
    distances=None,
    bidirectional=True,
    force2D=False,
    force2Ddimension=0,
    dtype=torch.int64,
    device=None,
):
    """
    Torch 版：行为与 cmatrices_generate_angles_np / C 实现对齐。

    参数：
        size: 1D tensor 或 array-like
        distances: 1D tensor 或 array-like 或 None
        bidirectional, force2D, force2Ddimension: 语义同上
        dtype: 输出 dtype（整数）
        device: 输出 device（如果 size / distances 是 tensor，会自动继承）
    """
    # 处理 device
    if torch.is_tensor(size):
        if device is None:
            device = size.device
        size_t = size.to(device=device, dtype=torch.int64).view(-1)
    else:
        device = device or "cpu"
        size_t = torch.as_tensor(size, device=device, dtype=torch.int64).view(-1)

    if size_t.ndim != 1:
        raise ValueError("Expected size to be 1D")
    Nd = int(size_t.numel())

    if distances is None:
        dist_t = torch.tensor([1], device=device, dtype=torch.int64)
    else:
        if torch.is_tensor(distances):
            dist_t = distances.to(device=device, dtype=torch.int64).view(-1)
        else:
            dist_t = torch.as_tensor(distances, device=device, dtype=torch.int64).view(-1)

    if dist_t.ndim != 1:
        raise ValueError("Expecting distances to be 1D")
    Ndist = int(dist_t.numel())

    # force2D 逻辑
    force2Ddim = int(force2Ddimension) if force2D else -1

    size_list = [int(x) for x in size_t.tolist()]
    dist_list = [int(x) for x in dist_t.tolist()]

    Na = _get_angle_count_core(size_list, dist_list, bidirectional=bool(bidirectional), force2Ddim=force2Ddim)
    if Na == 0:
        raise RuntimeError("Error getting angle count.")

    angles = _build_angles_core_to_tensor(
        size_list, dist_list, Nd, Ndist, force2Ddim, Na, device=device, dtype=dtype
    )

    return angles


# -----------------------------------------------------

import numpy as np


def get_angle_count_np(size, distances, bidirectional=True, force2Ddim=-1):
    """
    直接翻译 C 版 get_angle_count。
    size: 1D int array-like, shape [Nd]
    distances: 1D int array-like, shape [Ndist]
    """
    size = np.asarray(size, dtype=np.int64)
    distances = np.asarray(distances, dtype=np.int64)

    if size.ndim != 1:
        raise ValueError("Expected size to be 1D")
    if distances.ndim != 1:
        raise ValueError("Expected distances to be 1D")

    Nd = int(size.shape[0])
    Ndist = int(distances.shape[0])

    Na = 0
    for dist_idx in range(Ndist):
        d = int(distances[dist_idx])
        if d < 1:
            # C 里直接 return 0 表示错误
            return 0

        Na_d = 1
        Na_dd = 1
        for dim_idx in range(Nd):
            # 不在 out-of-plane 维度上生成角度
            if dim_idx == force2Ddim:
                continue

            if d < size[dim_idx]:
                # 全距离可用：±d
                Na_d *= (2 * d + 1)
                Na_dd *= (2 * d - 1)
            else:
                # 距离超过图像大小，只能用 size-1
                max_step = 2 * (int(size[dim_idx]) - 1) + 1
                Na_d *= max_step
                Na_dd *= max_step

        Na += (Na_d - Na_dd)

    if not bidirectional:
        Na //= 2

    return Na


def build_angles_np(size, distances, Nd, Ndist, force2Ddim, Na):
    """
    直接翻译 C 版 build_angles，返回 shape = (Na, Nd) 的 int 数组。

    注意：假定 size/distances/force2Ddim 与 get_angle_count_np 保持一致，
    且 Na 已经由 get_angle_count_np 计算好（包含 bidirectional 逻辑）。
    """
    size = np.asarray(size, dtype=np.int64)
    distances = np.asarray(distances, dtype=np.int64)

    size_list = [int(s) for s in size.tolist()]
    dist_list = [int(d) for d in distances.tolist()]

    # 找最大距离
    max_distance = 0
    for d in dist_list:
        if d < 1:
            # C 里返回 1 表示错误，这里直接抛异常更直观
            raise ValueError("Invalid distance (< 1) encountered in build_angles_np")
        if max_distance < d:
            max_distance = d

    n_offsets = 2 * max_distance + 1

    # offset_stride[Nd-1] = 1; 其余反向累乘
    offset_stride = np.empty(Nd, dtype=np.int64)
    offset_stride[Nd - 1] = 1
    for dim_idx in range(Nd - 2, -1, -1):
        offset_stride[dim_idx] = offset_stride[dim_idx + 1] * n_offsets

    # 预分配角度数组
    angles = np.empty((Na, Nd), dtype=np.int64)

    new_a_idx = 0  # 控制 offset 组合的计数器
    a_idx = 0      # 在 angles 中当前要填充的行

    while a_idx < Na:
        a_dist = 0       # 当前角度的 ∞ 范数
        valid_angle = True

        for dim_idx in range(Nd):
            # C: offset = max_distance - (new_a_idx / stride) % n_offsets;
            offset = max_distance - (
                (new_a_idx // int(offset_stride[dim_idx])) % n_offsets
            )

            # 条件同 C 版：
            # 1. force2D 维度上 offset 必须是 0
            # 2. offset 范围不能超出 [-(size[d]-1), (size[d]-1)]
            if ((dim_idx == force2Ddim and offset != 0) or
                offset >= size_list[dim_idx] or
                offset <= -size_list[dim_idx]):
                a_dist = -1  # 标记为非法
                valid_angle = False
                break

            # 写入当前候选角度（只有最终被确认“合法并需要”时才会推进 a_idx）
            angles[a_idx, dim_idx] = offset

            # 更新 ∞ 范数
            if a_dist < offset:
                a_dist = offset
            elif a_dist < -offset:
                a_dist = -offset

        new_a_idx += 1  # 不管角度合法与否，组合计数器都向前走

        # a_dist < 1: 非法 (-1) 或者是 (0,0,0)
        if a_dist < 1:
            continue

        if not valid_angle:
            continue

        # 检查此距离是否在 distances 中
        for d in dist_list:
            if a_dist == d:
                a_idx += 1  # 接受该角度
                break

    return angles


def generate_angles_np(
    size,
    distances=None,
    bidirectional=True,
    force2D=False,
    force2Ddimension=0,
):
    """
    Python/NumPy 版：重现 C 的 cmatrices_generate_angles 行为。

    参数：
        size: 1D array-like, 图像尺寸
        distances: 1D array-like 或 None。None 时默认 [1]
        bidirectional: 是否双向（True/False）
        force2D: 是否强制 2D
        force2Ddimension: 哪一维当作“out-of-plane”

    返回：
        numpy.ndarray, shape = (Na, Nd), dtype=int32（匹配 NPY_INT）
    """
    size_arr = np.asarray(size, dtype=np.int64)
    if size_arr.ndim != 1:
        raise ValueError("Expected a 1D array for size")
    Nd = int(size_arr.shape[0])

    if distances is None:
        distances_arr = np.array([1], dtype=np.int64)
    else:
        distances_arr = np.asarray(distances, dtype=np.int64)
        if distances_arr.ndim != 1:
            raise ValueError("Expecting distances array to be 1D")
    Ndist = int(distances_arr.shape[0])

    # C 里：如果不 force2D，则强制维度 = -1
    force2Ddim = int(force2Ddimension) if force2D else -1

    Na = get_angle_count_np(size_arr, distances_arr, bidirectional=bool(bidirectional), force2Ddim=force2Ddim)
    if Na == 0:
        raise RuntimeError("Error getting angle count.")

    angles = build_angles_np(size_arr, distances_arr, Nd, Ndist, force2Ddim, Na)

    # C 里是 NPY_INT，一般对应 C 的 int（常见是 int32）
    return angles.astype(np.int32, copy=False)

### glszm

def calculate_glszm_torch(
    image: torch.Tensor,
    mask: torch.Tensor,
    Ng: int,
    Ns: int,
    force2D: bool,
    force2Ddimension: int,
    kernelRadius: int = 0,
    voxels: torch.Tensor | None = None,
) -> torch.Tensor:
    """
    Torch 版本 GLSZM 计算核心，等价于 cmatrices_calculate_glszm + calculate_glszm + fill_glszm。

    参数
    ----
    image: Tensor, shape (Z,Y,X) 或 (Y,X)，已做 binning，类型任意（会转 int32）
    mask : Tensor, 同 shape，>0 视为在 ROI 内
    Ng   : 灰度级数量（最大 bin index）
    Ns   : ROI 内体素数（或 Nkernel 上限之前的 Ns），作为 zone-size 维度最大长度上界
    force2D, force2Ddimension, kernelRadius, voxels:
           同原 C 接口语义

    返回
    ----
    P_glszm: Tensor[float64], shape (Nvox, Ng, S)
             这里 S = Ns（上界），后续会在 RadiomicsGLSZM 里删掉全 0 列，得到真实 zone-size 轴长度。
    """
    device = image.device
    # 保证基本 dtype / 连续性
    image = image.to(torch.int32).contiguous()
    # mask 工作副本（会被本函数修改，外部不受影响）
    mask_work = (mask > 0).to(torch.bool).contiguous()

    Nd = image.dim()
    assert Nd in (2, 3), f"GLSZM 目前只考虑 2D/3D, got Nd={Nd}"
    size = list(image.shape)  # [Z,Y,X] or [Y,X]

    # ---- Nvox / voxel-based 相关 ----
    if voxels is not None:
        voxels = voxels.to(torch.long).contiguous()  # (Nd, Nvox)
        Nvox = int(voxels.shape[1])
    else:
        voxels = None
        Nvox = 1

    # 如果不强制 2D，则和 C 里一样把 force2Ddimension 置为 -1
    if not force2D:
        force2Ddimension = -1

    # ---- 计算 Nkernel 并调整 Ns（完全照抄 C 逻辑）----
    Ns_eff = Ns
    if voxels is not None and kernelRadius > 0:
        # v = Nd or Nd-1（如果 force2D）
        v_dim = Nd - (1 if force2D else 0)
        kernel_side = kernelRadius * 2 + 1
        Nkernel = kernel_side ** v_dim
        if Ns_eff > Nkernel:
            Ns_eff = Nkernel
    # 我们直接用 Ns_eff 作为 zone-size 维度上界
    max_region_dim = int(Ns_eff if Ns_eff > 0 else 1)

    # ---- 角度（邻居 offset）: build_angles_arr(NULL, ..., bidirectional=1）----
    angles, Na = build_angles_arr_torch(
        distances_obj=None,
        size=size,
        Nd=Nd,
        force2Ddimension=force2Ddimension,
        bidirectional=True,   # C 里最后一个参数 = 1
    )
    angles = angles.to(torch.int32).contiguous()  # (Na, Nd)

    # ---- 输出 GLSZM: (Nvox, Ng, max_region_dim) ----
    P_glszm = torch.zeros(
        (Nvox, Ng, max_region_dim),
        dtype=torch.float64,
        device=device,
    )

    # ---- 生成邻居 offset 列表（直接用 angles）----
    # angles[a, d] 就是第 a 个方向在第 d 维的偏移
    offsets = angles.tolist()  # List[List[int]]，Na x Nd

    # ---- region growing 主循环 ----
    # 如果是 voxel-based，要在每个 kernel 结束后恢复 mask_work 中的 1（与 C 里的 processedStack 逻辑一致）
    for v in range(Nvox):
        # ---- bounding box: set_bb_torch 和 C 的 set_bb 等价 ----
        bb_lo, bb_hi = set_bb_torch(
            v=v,
            size=size,
            voxels=voxels,
            Nd=Nd,
            Nvox=Nvox,
            kernelRadius=kernelRadius,
            force2Ddimension=force2Ddimension,
        )
        # processed 列表（只在 Nvox>1 时需要；Nvox=1 时不恢复 mask）
        processed_indices = []  # 存 (z,y,x) 或 (y,x)

        # 遍历 bounding box 中的每个 voxel，寻找“未处理 + 在 mask 中”的起始 voxel
        if Nd == 3:
            z0, y0, x0 = bb_lo
            z1, y1, x1 = bb_hi
            for z in range(z0, z1 + 1):
                for y in range(y0, y1 + 1):
                    for x in range(x0, x1 + 1):
                        if not mask_work[z, y, x]:
                            continue

                        # 当前区域起点
                        gl = int(image[z, y, x].item())
                        if gl <= 0 or gl > Ng:
                            # 非法灰度直接跳过（C 里写 gl<=0 会在 fill_glszm 时触发报错，这里直接不计）
                            mask_work[z, y, x] = False
                            if Nvox > 1:
                                processed_indices.append((z, y, x))
                            continue

                        region_size = 0
                        # BFS / DFS 的栈
                        stack = [(z, y, x)]

                        # 起点标记为已处理
                        mask_work[z, y, x] = False
                        if Nvox > 1:
                            processed_indices.append((z, y, x))

                        while stack:
                            cz, cy, cx = stack.pop()
                            region_size += 1

                            # 遍历邻居
                            for off in offsets:
                                dz, dy, dx = off
                                nz, ny, nx = cz + dz, cy + dy, cx + dx

                                # 在 bounding box 内？
                                if (
                                    nz < z0 or nz > z1 or
                                    ny < y0 or ny > y1 or
                                    nx < x0 or nx > x1
                                ):
                                    continue

                                if not mask_work[nz, ny, nx]:
                                    continue

                                if int(image[nz, ny, nx].item()) != gl:
                                    continue

                                mask_work[nz, ny, nx] = False
                                stack.append((nz, ny, nx))
                                if Nvox > 1:
                                    processed_indices.append((nz, ny, nx))

                        # 一个 region 完成，填入 GLSZM
                        if region_size <= 0:
                            continue
                        # C 里 maxSize <= Ns_eff，理论上 region_size 不会大于 Ns_eff
                        if region_size > max_region_dim:
                            # 理论上不会走到这里；为安全起见做一个 clamp
                            region_size = max_region_dim

                        # 累加：gray level gl → index gl-1; size j → index j-1
                        P_glszm[v, gl - 1, region_size - 1] += 1

        elif Nd == 2:
            y0, x0 = bb_lo
            y1, x1 = bb_hi
            for y in range(y0, y1 + 1):
                for x in range(x0, x1 + 1):
                    if not mask_work[y, x]:
                        continue

                    gl = int(image[y, x].item())
                    if gl <= 0 or gl > Ng:
                        mask_work[y, x] = False
                        if Nvox > 1:
                            processed_indices.append((y, x))
                        continue

                    region_size = 0
                    stack = [(y, x)]

                    mask_work[y, x] = False
                    if Nvox > 1:
                        processed_indices.append((y, x))

                    while stack:
                        cy, cx = stack.pop()
                        region_size += 1

                        for off in offsets:
                            dy, dx = off
                            ny, nx = cy + dy, cx + dx
                            if (
                                ny < y0 or ny > y1 or
                                nx < x0 or nx > x1
                            ):
                                continue

                            if not mask_work[ny, nx]:
                                continue

                            if int(image[ny, nx].item()) != gl:
                                continue

                            mask_work[ny, nx] = False
                            stack.append((ny, nx))
                            if Nvox > 1:
                                processed_indices.append((ny, nx))

                    if region_size <= 0:
                        continue
                    if region_size > max_region_dim:
                        region_size = max_region_dim

                    P_glszm[v, gl - 1, region_size - 1] += 1

        # voxel-based: 恢复 mask（完全对应 C 里的 processedStack 回滚）
        if Nvox > 1 and processed_indices:
            if Nd == 3:
                for (z, y, x) in processed_indices:
                    mask_work[z, y, x] = True
            else:
                for (y, x) in processed_indices:
                    mask_work[y, x] = True

    # 返回形状 (Nvox, Ng, max_region_dim)
    return P_glszm

def calculate_glszm_torch_vectorized(
    image: torch.Tensor,
    mask: torch.Tensor,
    Ng: int,
    Ns: int,
    force2D: bool,
    force2Ddimension: int,
    kernelRadius: int = 0,
    voxels: torch.Tensor | None = None,
) -> torch.Tensor:
    """
    向量化 / 卷积版 GLSZM 计算（ROI 级别），对应原始 C 版本：
      cmatrices_calculate_glszm + calculate_glszm + fill_glszm
    但只实现了非 voxel-based 的整块 ROI 情况 (Nvox = 1)。

    参数
    ----
    image: Tensor, shape (Z, Y, X) 或 (Y, X)，已经 binning，整型值 1..Ng
    mask : Tensor, 同 shape，>0 的位置视为在 ROI 内
    Ng   : 灰度级数量（最大离散 bin）
    Ns   : ROI 体素数（或者你在 Python 端传入的 Ns，上界，用来做 zone-size 轴长度）
    force2D: 目前要求 False（3D 用 26-邻域，2D 用 8-邻域）
    force2Ddimension: 暂时忽略（因为 force2D=False）
    kernelRadius, voxels: 目前不支持 voxel-based，要求 kernelRadius=0, voxels=None

    返回
    ----
    P_glszm: Tensor[float64], shape (1, Ng, Ns)
             后续你会在 RadiomicsGLSZM 里删除全 0 的 size 列，得到真实 Ns。
    """
    if voxels is not None or kernelRadius != 0:
        raise NotImplementedError(
            "calculate_glszm_torch_vectorized 目前只支持整块 ROI (voxels=None, kernelRadius=0)"
        )
    if force2D:
        raise NotImplementedError(
            "向量化版本暂未实现 force2D=True 的行为（暂支持 3D/2D 全 26/8 连通）"
        )

    device = image.device

    # 保证 dtype / 连续性
    image = image.to(torch.int64).contiguous()
    mask_work = (mask > 0).to(torch.bool).contiguous()

    Nd = image.dim()
    assert Nd in (2, 3), f"GLSZM 目前只考虑 2D/3D, got Nd={Nd}"

    # 输出 GLSZM: Nvox=1（整块 ROI）
    P_glszm = torch.zeros((1, Ng, Ns), dtype=torch.float64, device=device)

    # 卷积核（3x3 或 3x3x3，全 1，对应 Chebyshev 距离 1 的邻域）
    if Nd == 3:
        kernel = torch.ones((1, 1, 3, 3, 3), dtype=torch.float32, device=device)
    else:  # Nd == 2
        kernel = torch.ones((1, 1, 3, 3), dtype=torch.float32, device=device)

    # ======================
    # 主循环：不断从 ROI 中拿一个种子，做连通区域扩展
    # ======================
    # 注意：这里没有 per-voxel 循环，而是：
    #  - 每次用 conv2d/conv3d 做 region dilation
    #  - 只在 Python 层按「区域个数」循环
    while mask_work.any():
        # 选一个未处理的 seed voxel
        seed_idx = mask_work.nonzero(as_tuple=False)[0]  # (Nd,)
        if Nd == 3:
            z0, y0, x0 = seed_idx.tolist()
            gl = int(image[z0, y0, x0].item())
        else:
            y0, x0 = seed_idx.tolist()
            gl = int(image[y0, x0].item())

        # 非法灰度直接标记处理掉
        if gl <= 0 or gl > Ng:
            if Nd == 3:
                mask_work[z0, y0, x0] = False
            else:
                mask_work[y0, x0] = False
            continue

        # 初始化当前 region 掩码（仅 seed = True）
        region = torch.zeros_like(mask_work)
        if Nd == 3:
            region[z0, y0, x0] = True
        else:
            region[y0, x0] = True

        # ----------------------
        # 用 conv 反复扩张 region，直到不再增长
        # ----------------------
        while True:
            # [B=1,C=1,...] 形式做卷积
            region_f = region.to(torch.float32).unsqueeze(0).unsqueeze(0)  # (1,1,*,*,*)
            if Nd == 3:
                dilated = F.conv3d(region_f, kernel, padding=1)
            else:
                dilated = F.conv2d(region_f, kernel, padding=1)

            dilated = dilated.squeeze(0).squeeze(0) > 0

            # 只能在 ROI 且 和 seed 同灰度，且未被其他 region 占用的位置扩张
            same_gray = (image == gl)
            grown = region | (dilated & mask_work & same_gray)

            if torch.equal(grown, region):
                break
            region = grown

        # 得到一个完整 zone，计算大小
        size_zone = int(region.sum().item())
        if size_zone <= 0:
            # 理论上不会发生，但安全起见
            if Nd == 3:
                mask_work[z0, y0, x0] = False
            else:
                mask_work[y0, x0] = False
            continue

        if size_zone > Ns:
            # 与 C 版本一致：理论上不会 > Ns，这里 clamp 一下防炸
            size_zone = Ns

        # 记录到 GLSZM: gl -> index gl-1, size_zone -> index size_zone-1
        P_glszm[0, gl - 1, size_zone - 1] += 1

        # 当前 region 所有 voxel 标记为已处理（不再作为后续种子）
        mask_work[region] = False

    return P_glszm

### ngtdm

# cmatrices.py
# def calculate_ngtdm_torch(
#     image: torch.Tensor,
#     mask: torch.Tensor,
#     distances,
#     Ng: int,
#     force2D: bool,
#     force2Ddimension: int,
#     kernelRadius: int = 0,
#     voxels: torch.Tensor | None = None,
#     device=None,
#     dtype=torch.float32,
# ) -> torch.Tensor:
#     """
#     Torch 版 NGTDM 计算，逻辑对应 C 里的 cmatrices_calculate_ngtdm + calculate_ngtdm。

#     参数
#     ----
#     image : torch.Tensor (int)
#         N 维整型张量，形状 (N1, N2, ..., Nd)，灰度值范围 [1..Ng]。
#     mask : torch.Tensor (bool 或 0/1)
#         同 shape，表示 ROI。
#     distances : 1D list/tensor
#         距离列表，对应 C 版 distances_obj。
#     Ng : int
#         灰度级数。
#     force2D : bool
#         是否强制 2D。
#     force2Ddimension : int
#         如果 force2D=True，则这个维度视为 out-of-plane（只允许 offset=0）。
#     kernelRadius : int, default 0
#         体素级（voxel-based）提取的核半径；0 表示整块 ROI（非 voxel-based）。
#     voxels : torch.Tensor | None
#         体素级提取时的中心点索引，形状 (Nvox, Nd)，每行是一个体素的坐标。
#     device : torch.device | None
#         计算使用的 device，默认 image.device。
#     dtype : torch.dtype
#         输出的浮点精度，默认 torch.double。

#     返回
#     ----
#     torch.Tensor
#         形状 (Nvox, Ng, 3)，其中最后一维分别为:
#         [:, :, 0] = n_i
#         [:, :, 1] = s_i
#         [:, :, 2] = 灰度值 i (1..Ng)
#     """
#     if device is None:
#         device = image.device

#     image = image.to(device=device, dtype=torch.long)
#     mask = mask.to(device=device)
#     if mask.dtype != torch.bool:
#         mask = mask != 0

#     Nd = image.ndim
#     size = [int(s) for s in image.shape]

#     # 处理 distances
#     distances = torch.as_tensor(distances, device=device, dtype=torch.long).view(-1)
#     distances_list = [int(d.item()) for d in distances]
#     if len(distances_list) == 0:
#         raise ValueError("distances must contain at least one element")

#     # force2D 逻辑：和 C 里一样，不 force2D 就把维度设为 -1
#     if not force2D:
#         force2Ddimension = -1
#     else:
#         if not (0 <= force2Ddimension < Nd):
#             force2Ddimension = -1

#     # ---- 下面两个 helper 是 C 版 get_angle_count / build_angles 的 torch 版本 ----

#     def _get_angle_count(size_, distances_, Nd_, bidirectional: bool, force2Ddim: int) -> int:
#         Na = 0
#         for dist in distances_:
#             if dist < 1:
#                 return 0
#             Na_d = 1
#             Na_dd = 1
#             for dim_idx in range(Nd_):
#                 if dim_idx == force2Ddim:
#                     continue
#                 if dist < size_[dim_idx]:
#                     Na_d *= (2 * dist + 1)
#                     Na_dd *= (2 * dist - 1)
#                 else:
#                     # 被图像尺寸限制
#                     Na_d *= (2 * (size_[dim_idx] - 1) + 1)
#                     Na_dd *= (2 * (size_[dim_idx] - 1) + 1)
#             Na += (Na_d - Na_dd)
#         # NGTDM 用的是 bidirectional=True，所以不会走到这里
#         if not bidirectional:
#             Na //= 2
#         return Na

#     def _build_angles(size_, distances_, Nd_, force2Ddim: int, bidirectional: bool = True) -> torch.Tensor:
#         size_ = list(size_)
#         distances_ = list(distances_)
#         Na = _get_angle_count(size_, distances_, Nd_, bidirectional, force2Ddim)
#         if Na <= 0:
#             raise ValueError("No valid angles could be generated for the given distances and image size")

#         max_distance = max(distances_)
#         n_offsets = 2 * max_distance + 1

#         # offset_stride 用来枚举各维 offset 组合
#         offset_stride = [0] * Nd_
#         offset_stride[Nd_ - 1] = 1
#         for dim_idx in range(Nd_ - 2, -1, -1):
#             offset_stride[dim_idx] = offset_stride[dim_idx + 1] * n_offsets

#         angles = [[0] * Nd_ for _ in range(Na)]
#         new_a_idx = 0  # 控制 offset 组合
#         a_idx = 0      # 已接受的 angle 数

#         while a_idx < Na:
#             a_dist = 0
#             candidate = [0] * Nd_
#             for dim_idx in range(Nd_):
#                 offset = max_distance - (new_a_idx // offset_stride[dim_idx]) % n_offsets
#                 # 非法 angle：超出尺寸或在 force2Ddim 上有非 0 offset
#                 if (
#                     (dim_idx == force2Ddim and offset != 0)
#                     or offset >= size_[dim_idx]
#                     or offset <= -size_[dim_idx]
#                 ):
#                     a_dist = -1
#                     break
#                 candidate[dim_idx] = offset
#                 if a_dist < abs(offset):
#                     a_dist = abs(offset)

#             new_a_idx += 1

#             # a_dist < 1: 要么 0 向量，要么非法 angle，直接丢弃
#             if a_dist < 1:
#                 continue

#             # 只保留距离在 distances_ 里的 angle
#             if any(a_dist == d for d in distances_):
#                 angles[a_idx] = candidate
#                 a_idx += 1

#         return torch.tensor(angles, dtype=torch.long, device=device)

#     # 为 NGTDM 生成 angles（bidirectional=True）
#     angles = _build_angles(size, distances_list, Nd, force2Ddimension, bidirectional=True)
#     Na = int(angles.shape[0])

#     # 按 C 顺序计算展平 strides（元素数量，而不是字节）
#     strides = [1] * Nd
#     for d in range(Nd - 2, -1, -1):
#         strides[d] = strides[d + 1] * size[d + 1]

#     image_flat = image.reshape(-1)
#     mask_flat = mask.reshape(-1)

#     def _single_ngtdm(bb: list[int]) -> torch.Tensor:
#         """
#         对一个 bounding box 计算 NGTDM。
#         bb: 长度 2*Nd 的列表 [lo0..lo{Nd-1}, hi0..hi{Nd-1}]
#         """
#         ngtdm = torch.zeros((Ng, 3), dtype=dtype, device=device)
#         # 第 2 列填 1..Ng
#         ngtdm[:, 2] = torch.arange(1, Ng + 1, dtype=dtype, device=device)

#         # 整个 image 的元素个数
#         Ni = 1
#         for s in size:
#             Ni *= s

#         # 起始 flat index（bb 的下界）
#         i = 0
#         for d in range(Nd):
#             i += bb[d] * strides[d]

#         cur_idx = [0] * Nd
#         ngtdm_idx_max = Ng * 3

#         while i < Ni:
#             # 先在各维上把 i 调到 bounding box 内
#             for d in range(Nd - 1, 0, -1):
#                 cur_idx[d] = (i % strides[d - 1]) // strides[d]
#                 if cur_idx[d] > bb[Nd + d]:
#                     i += (size[d] - cur_idx[d] + bb[d]) * strides[d]
#                     cur_idx[d] = bb[d]
#                 elif cur_idx[d] < bb[d]:
#                     i += (bb[d] - cur_idx[d]) * strides[d]
#                     cur_idx[d] = bb[d]

#             cur_idx[0] = i // strides[0]
#             if cur_idx[0] > bb[Nd]:
#                 # 超出第 0 维上界，整个 bounding box 结束
#                 break

#             if mask_flat[i]:
#                 count = 0.0
#                 ssum = 0.0

#                 # 遍历所有邻域方向
#                 for a in range(Na):
#                     j = i
#                     for d in range(Nd):
#                         off = int(angles[a, d].item())
#                         idx_d = cur_idx[d]
#                         # 出界就标记 j=i，表示该方向无邻居
#                         if idx_d + off < bb[d] or idx_d + off > bb[Nd + d]:
#                             j = i
#                             break
#                         j += off * strides[d]

#                     if j != i and mask_flat[j]:
#                         count += 1.0
#                         ssum += int(image_flat[j].item())

#                 if count == 0.0:
#                     diff = 0.0
#                 else:
#                     diff = float(int(image_flat[i].item()) - (ssum / count))
#                     if diff < 0.0:
#                         diff = -diff

#                 gl_val = int(image_flat[i].item())
#                 ngtdm_idx = (gl_val - 1) * 3
#                 if gl_val <= 0 or ngtdm_idx >= ngtdm_idx_max:
#                     raise IndexError("NGTDM: gray level index out of range")

#                 # n_i
#                 ngtdm[gl_val - 1, 0] += 1.0
#                 # s_i 累加绝对差
#                 ngtdm[gl_val - 1, 1] += diff

#             i += 1

#         return ngtdm

#     # ---- voxel-based / 非 voxel-based 两种模式 ----

#     if voxels is None or kernelRadius <= 0:
#         # 整块 ROI，一个 NGTDM
#         bb_lo = [0] * Nd
#         bb_hi = [s - 1 for s in size]
#         bb = bb_lo + bb_hi
#         result = _single_ngtdm(bb).unsqueeze(0)  # (1, Ng, 3)
#     else:
#         voxels = voxels.to(device=device, dtype=torch.long)
#         if voxels.ndim != 2 or voxels.shape[1] != Nd:
#             raise ValueError(f"voxels must have shape (Nvox, Nd={Nd}), got {tuple(voxels.shape)}")

#         Nvox = int(voxels.shape[0])
#         result = torch.empty((Nvox, Ng, 3), dtype=dtype, device=device)

#         for v in range(Nvox):
#             bb_lo: list[int] = []
#             bb_hi: list[int] = []
#             for d in range(Nd):
#                 center_d = int(voxels[v, d].item())
#                 if force2D and d == force2Ddimension:
#                     lo = hi = center_d
#                 else:
#                     lo = max(center_d - kernelRadius, 0)
#                     hi = min(center_d + kernelRadius, size[d] - 1)
#                 bb_lo.append(lo)
#                 bb_hi.append(hi)
#             bb = bb_lo + bb_hi
#             result[v] = _single_ngtdm(bb)

#     return result

def calculate_ngtdm_torch(
    image: torch.Tensor,
    mask: torch.Tensor,
    distances,
    Ng: int,
    force2D: bool,
    force2Ddimension: int,
    kernelRadius: int = 0,
    voxels: Optional[torch.Tensor] = None,
    device=None,
    dtype: torch.dtype = torch.float32,
) -> torch.Tensor:
    """
    Torch 版 NGTDM 计算，逻辑对应 C 里的 cmatrices_calculate_ngtdm + calculate_ngtdm。

    参数
    ----
    image : torch.Tensor (int)
        N 维整型张量，形状 (N1, N2, ..., Nd)，灰度值范围 [1..Ng]。
    mask : torch.Tensor (bool 或 0/1)
        同 shape，表示 ROI。
    distances : 1D list/tensor
        距离列表，对应 C 版 distances_obj。
    Ng : int
        灰度级数。
    force2D : bool
        是否强制 2D。
    force2Ddimension : int
        如果 force2D=True，则这个维度视为 out-of-plane（只允许 offset=0）。
    kernelRadius : int, default 0
        体素级（voxel-based）提取的核半径；0 表示整块 ROI（非 voxel-based）。
    voxels : torch.Tensor | None
        体素级提取时的中心点索引，形状 (Nvox, Nd)，每行是一个体素的坐标。
    device : torch.device | None
        计算使用的 device，默认 image.device。
    dtype : torch.dtype
        输出的浮点精度，默认 torch.float32。

    返回
    ----
    torch.Tensor
        形状 (Nvox, Ng, 3)，其中最后一维分别为:
        [:, :, 0] = n_i
        [:, :, 1] = s_i
        [:, :, 2] = 灰度值 i (1..Ng)
    """
    if device is None:
        device = image.device

    image = image.to(device=device, dtype=torch.long)
    mask = mask.to(device=device)
    if mask.dtype != torch.bool:
        mask = mask != 0

    Nd = image.ndim
    size: List[int] = [int(s) for s in image.shape]

    # 处理 distances
    distances = torch.as_tensor(distances, device=device, dtype=torch.long).view(-1)
    distances_list: List[int] = [int(d.item()) for d in distances]
    if len(distances_list) == 0:
        raise ValueError("distances must contain at least one element")

    # force2D 逻辑：和 C 一样，不 force2D 就把维度设为 -1
    if not force2D:
        force2Ddimension = -1
    else:
        if not (0 <= force2Ddimension < Nd):
            force2Ddimension = -1

    # ================== helper: angle count & build angles ==================

    def _get_angle_count(
        size_: List[int],
        distances_: List[int],
        Nd_: int,
        bidirectional: bool,
        force2Ddim: int,
    ) -> int:
        Na = 0
        for dist in distances_:
            if dist < 1:
                return 0
            Na_d = 1
            Na_dd = 1
            for dim_idx in range(Nd_):
                if dim_idx == force2Ddim:
                    continue
                if dist < size_[dim_idx]:
                    Na_d *= (2 * dist + 1)
                    Na_dd *= (2 * dist - 1)
                else:
                    # 被图像尺寸限制
                    Na_d *= (2 * (size_[dim_idx] - 1) + 1)
                    Na_dd *= (2 * (size_[dim_idx] - 1) + 1)
            Na += (Na_d - Na_dd)
        # NGTDM 用的是 bidirectional=True，所以不会走到这里
        if not bidirectional:
            Na //= 2
        return Na

    def _build_angles(
        size_: List[int],
        distances_: List[int],
        Nd_: int,
        force2Ddim: int,
        bidirectional: bool = True,
    ) -> torch.Tensor:
        size_ = list(size_)
        distances_ = list(distances_)
        Na = _get_angle_count(size_, distances_, Nd_, bidirectional, force2Ddim)
        if Na <= 0:
            raise ValueError("No valid angles could be generated for the given distances and image size")

        max_distance = max(distances_)
        n_offsets = 2 * max_distance + 1

        # offset_stride 用来枚举各维 offset 组合
        offset_stride = [0] * Nd_
        offset_stride[Nd_ - 1] = 1
        for dim_idx in range(Nd_ - 2, -1, -1):
            offset_stride[dim_idx] = offset_stride[dim_idx + 1] * n_offsets

        angles = [[0] * Nd_ for _ in range(Na)]
        new_a_idx = 0  # 控制 offset 组合
        a_idx = 0      # 已接受的 angle 数

        while a_idx < Na:
            a_dist = 0
            candidate = [0] * Nd_
            for dim_idx in range(Nd_):
                offset = max_distance - (new_a_idx // offset_stride[dim_idx]) % n_offsets
                # 非法 angle：超出尺寸或在 force2Ddim 上有非 0 offset
                if (
                    (dim_idx == force2Ddim and offset != 0)
                    or offset >= size_[dim_idx]
                    or offset <= -size_[dim_idx]
                ):
                    a_dist = -1
                    break
                candidate[dim_idx] = offset
                if a_dist < abs(offset):
                    a_dist = abs(offset)

            new_a_idx += 1

            # a_dist < 1: 要么 0 向量，要么非法 angle，直接丢弃
            if a_dist < 1:
                continue

            # 只保留距离在 distances_ 里的 angle
            if any(a_dist == d for d in distances_):
                angles[a_idx] = candidate
                a_idx += 1

        return torch.tensor(angles, dtype=torch.long, device=device)

    # 为 NGTDM 生成 angles（bidirectional=True）
    angles = _build_angles(size, distances_list, Nd, force2Ddimension, bidirectional=True)
    Na = int(angles.shape[0])

    # 按 C 顺序计算展平 strides（元素数量，而不是字节）
    strides: List[int] = [1] * Nd
    for d in range(Nd - 2, -1, -1):
        strides[d] = strides[d + 1] * size[d + 1]

    strides_t = torch.tensor(strides, device=device, dtype=torch.long)

    image_flat = image.reshape(-1)
    mask_flat = mask.reshape(-1)

    # 整个 image 元素个数（给 bbox 版本用）
    Ni_total = 1
    for s in size:
        Ni_total *= s

    # ================== 优化后的整块 ROI 版本（非 voxel-based） ==================
    def _single_ngtdm_vectorized() -> torch.Tensor:
        """
        整块 ROI 的 NGTDM（非 voxel-based），向量化实现。
        返回 shape (Ng, 3): [n_i, s_i, gray_level_i]
        """
        # 1. ROI 内所有体素坐标和线性索引
        coords = torch.nonzero(mask, as_tuple=False)  # (Nvox, Nd)
        if coords.numel() == 0:
            out = torch.zeros((Ng, 3), dtype=dtype, device=device)
            out[:, 2] = torch.arange(1, Ng + 1, dtype=dtype, device=device)
            return out

        coords = coords.to(torch.long)
        Nvox = coords.shape[0]

        base_idx = (coords * strides_t).sum(dim=1)  # (Nvox,)
        center_gray = image_flat[base_idx].to(torch.long)  # (Nvox,)

        count = torch.zeros(Nvox, dtype=dtype, device=device)
        ssum = torch.zeros(Nvox, dtype=dtype, device=device)

        size_t = torch.tensor(size, device=device, dtype=torch.long)

        # 2. 对每个 angle，向量化计算邻居并累积
        for a in range(Na):
            offset_vec = angles[a]  # (Nd,)
            # 邻居坐标
            neigh_coords = coords + offset_vec  # (Nvox, Nd)

            # 先做边界检查：0 <= coord < size
            valid = torch.ones(Nvox, dtype=torch.bool, device=device)
            for d in range(Nd):
                valid &= (neigh_coords[:, d] >= 0) & (neigh_coords[:, d] < size_t[d])

            if not valid.any():
                continue

            valid_idx = torch.nonzero(valid, as_tuple=False).squeeze(1)
            neigh_coords_valid = neigh_coords[valid_idx]
            neigh_idx_valid = (neigh_coords_valid * strides_t).sum(dim=1)

            # 邻居是否在 ROI 内
            neigh_in_roi = mask_flat[neigh_idx_valid]
            if not neigh_in_roi.any():
                continue

            neigh_valid_idx = valid_idx[neigh_in_roi]
            neigh_idx_roi = neigh_idx_valid[neigh_in_roi]

            # 有效邻居灰度值
            neigh_val = image_flat[neigh_idx_roi].to(dtype)

            # 对应中心体素累加
            count[neigh_valid_idx] += 1.0
            ssum[neigh_valid_idx] += neigh_val

        # 3. 计算每个体素的 diff = |center - neighbor_mean|
        diff = torch.zeros(Nvox, dtype=dtype, device=device)
        has_neighbors = count > 0
        if has_neighbors.any():
            neighbor_mean = torch.zeros_like(diff)
            neighbor_mean[has_neighbors] = ssum[has_neighbors] / count[has_neighbors]
            diff[has_neighbors] = (center_gray.to(dtype) - neighbor_mean[has_neighbors]).abs()

        # 4. 按灰度级聚合到 n_i & s_i
        gl = center_gray.clamp(min=1, max=Ng)
        gl_idx = gl - 1  # 0..Ng-1
        ones = torch.ones_like(diff)

        n_i = torch.bincount(gl_idx, weights=ones, minlength=Ng).to(dtype)
        s_i = torch.bincount(gl_idx, weights=diff, minlength=Ng).to(dtype)

        ngtdm = torch.empty((Ng, 3), dtype=dtype, device=device)
        ngtdm[:, 0] = n_i
        ngtdm[:, 1] = s_i
        ngtdm[:, 2] = torch.arange(1, Ng + 1, dtype=dtype, device=device)
        return ngtdm

    # ================== 原始 bbox 扫描版（给 voxel-based 用） ==================
    def _single_ngtdm_bbox(bb: List[int]) -> torch.Tensor:
        """
        对一个 bounding box 计算 NGTDM。
        bb: 长度 2*Nd 的列表 [lo0..lo{Nd-1}, hi0..hi{Nd-1}]
        使用的是你原来的 while + cur_idx + strides 扫描逻辑。
        """
        ngtdm = torch.zeros((Ng, 3), dtype=dtype, device=device)
        # 第 2 列填 1..Ng
        ngtdm[:, 2] = torch.arange(1, Ng + 1, dtype=dtype, device=device)

        ngtdm_idx_max = Ng * 3

        i = 0  # flat index 起点
        for d in range(Nd):
            i += bb[d] * strides[d]

        cur_idx = [0] * Nd

        while i < Ni_total:
            # 先在各维上把 i 调到 bounding box 内
            for d in range(Nd - 1, 0, -1):
                cur_idx[d] = (i % strides[d - 1]) // strides[d]
                if cur_idx[d] > bb[Nd + d]:
                    i += (size[d] - cur_idx[d] + bb[d]) * strides[d]
                    cur_idx[d] = bb[d]
                elif cur_idx[d] < bb[d]:
                    i += (bb[d] - cur_idx[d]) * strides[d]
                    cur_idx[d] = bb[d]

            cur_idx[0] = i // strides[0]
            if cur_idx[0] > bb[Nd]:
                # 超出第 0 维上界，整个 bounding box 结束
                break

            if mask_flat[i]:
                count = 0.0
                ssum = 0.0

                # 遍历所有邻域方向
                for a in range(Na):
                    j = i
                    for d in range(Nd):
                        off = int(angles[a, d].item())
                        idx_d = cur_idx[d]
                        # 出界就标记 j=i，表示该方向无邻居
                        if idx_d + off < bb[d] or idx_d + off > bb[Nd + d]:
                            j = i
                            break
                        j += off * strides[d]

                    if j != i and mask_flat[j]:
                        count += 1.0
                        ssum += int(image_flat[j].item())

                if count == 0.0:
                    diff = 0.0
                else:
                    diff = float(int(image_flat[i].item()) - (ssum / count))
                    if diff < 0.0:
                        diff = -diff

                gl_val = int(image_flat[i].item())
                ngtdm_idx = (gl_val - 1) * 3
                if gl_val <= 0 or ngtdm_idx >= ngtdm_idx_max:
                    raise IndexError("NGTDM: gray level index out of range")

                # n_i
                ngtdm[gl_val - 1, 0] += 1.0
                # s_i 累加绝对差
                ngtdm[gl_val - 1, 1] += diff

            i += 1

        return ngtdm

    # ================== 分支：非 voxel / voxel ==================
    if voxels is None or kernelRadius <= 0:
        # 整块 ROI，一个 NGTDM（优化后）
        result = _single_ngtdm_vectorized().unsqueeze(0)  # (1, Ng, 3)
    else:
        # voxel-based：逐 voxel 构造 bbox，使用原始 bbox 版本
        voxels = voxels.to(device=device, dtype=torch.long)
        if voxels.ndim != 2 or voxels.shape[1] != Nd:
            raise ValueError(f"voxels must have shape (Nvox, Nd={Nd}), got {tuple(voxels.shape)}")

        Nvox = int(voxels.shape[0])
        result = torch.empty((Nvox, Ng, 3), dtype=dtype, device=device)

        for v in range(Nvox):
            bb_lo: List[int] = []
            bb_hi: List[int] = []
            for d in range(Nd):
                center_d = int(voxels[v, d].item())
                if force2D and d == force2Ddimension:
                    lo = hi = center_d
                else:
                    lo = max(center_d - kernelRadius, 0)
                    hi = min(center_d + kernelRadius, size[d] - 1)
                bb_lo.append(lo)
                bb_hi.append(hi)
            bb = bb_lo + bb_hi
            result[v] = _single_ngtdm_bbox(bb)

    return result