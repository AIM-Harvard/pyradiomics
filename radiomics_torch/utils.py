import torch

def torch_delete(arr, obj, axis=None):
    """Delete sub-arrays from a tensor along a given axis.

    This function mimics the behavior of numpy.delete for PyTorch tensors.

    Parameters
    ----------
    arr : torch.Tensor
        Input tensor.
    obj : int, slice, or sequence of ints
        Indices of sub-arrays to delete.
    axis : int, optional
        Axis along which to delete the sub-arrays. If None, the input tensor is flattened.

    Returns
    -------
    torch.Tensor
        Tensor with specified sub-arrays deleted.
    """
    if axis is None:
        arr = arr.flatten()
        axis = 0

    mask = torch.ones(arr.size(axis), dtype=torch.bool, device=arr.device)
    if isinstance(obj, int):
        obj = [obj]
    mask[obj] = False

    return torch.index_select(arr, axis, torch.nonzero(mask, as_tuple=False).squeeze())

# ========================================
# NaN-safe ops with dim support (numpy-like)
# ========================================

def nanmax(tensor, dim=None, keepdim=False):
    min_value = torch.finfo(tensor.dtype).min
    output = tensor.nan_to_num(min_value).max(dim=dim, keepdim=keepdim)
    return output


def nanmin(tensor, dim=None, keepdim=False):
    max_value = torch.finfo(tensor.dtype).max
    output = tensor.nan_to_num(max_value).min(dim=dim, keepdim=keepdim)
    return output


def nanvar(tensor, dim=None, keepdim=False):
    tensor_mean = tensor.nanmean(dim=dim, keepdim=True)
    output = (tensor - tensor_mean).square().nanmean(dim=dim, keepdim=keepdim)
    return output


def nanstd(tensor, dim=None, keepdim=False):
    output = nanvar(tensor, dim=dim, keepdim=keepdim)
    output = output.sqrt()
    return output


def nanprod(tensor, dim=None, keepdim=False):
    output = tensor.nan_to_num(1).prod(dim=dim, keepdim=keepdim)
    return output


def nancumprod(tensor, dim=None, keepdim=False):
    output = tensor.nan_to_num(1).cumprod(dim=dim, keepdim=keepdim)
    return output


def nancumsum(tensor, dim=None, keepdim=False):
    output = tensor.nan_to_num(0).cumsum(dim=dim, keepdim=keepdim)
    return output


def nanargmin(tensor, dim=None, keepdim=False):
    max_value = torch.finfo(tensor.dtype).max
    output = tensor.nan_to_num(max_value).argmin(dim=dim, keepdim=keepdim)
    return output


def nanargmax(tensor, dim=None, keepdim=False):
    min_value = torch.finfo(tensor.dtype).min
    output = tensor.nan_to_num(min_value).argmax(dim=dim, keepdim=keepdim)
    return output

def nanquantile(tensor: torch.Tensor,
                q: float,
                dim: int | None = None,
                keepdim: bool = False) -> torch.Tensor:
    """
    PyTorch 版 nanquantile，行为类似 np.nanquantile
    - q: in [0, 1]，例如 0.1 / 0.9 / 0.25 / 0.75
    - dim=None: 在整个 tensor 上算一个标量
    - dim=int: 在指定维度上算分位数
    """
    if dim is None:
        # 全局把 NaN 去掉
        x = tensor[~torch.isnan(tensor)]
        if x.numel() == 0:
            return torch.tensor(float("nan"), device=tensor.device, dtype=tensor.dtype)
        return torch.quantile(x, q)

    # 统一正索引
    dim = dim if dim >= 0 else tensor.dim() + dim

    # 把要计算的 dim 挪到最后，方便展平
    x = tensor.transpose(dim, -1)        # shape: (..., K)
    K = x.shape[-1]
    x_flat = x.reshape(-1, K)            # (M, K)，每一行是一段要算 quantile 的数据

    res_list = []
    for row in x_flat:
        v = row[~torch.isnan(row)]
        if v.numel() == 0:
            res_list.append(torch.tensor(float("nan"), device=tensor.device,
                                         dtype=tensor.dtype))
        else:
            res_list.append(torch.quantile(v, q))
    res = torch.stack(res_list)          # (M,)

    # 还原成原来的 shape（在 dim 上 reduce 之后）
    orig_shape = list(tensor.shape)
    if keepdim:
        orig_shape[dim] = 1
    else:
        del orig_shape[dim]

    res = res.reshape(orig_shape)
    return res

def delete_torch(x: torch.Tensor, indices: torch.Tensor, dim: int = 1):
    indices = indices.to(dtype=torch.long, device=x.device)
    mask = torch.ones(x.size(dim), dtype=torch.bool, device=x.device)

    mask[indices] = False
    keep = mask.nonzero(as_tuple=True)[0]

    return x.index_select(dim, keep)