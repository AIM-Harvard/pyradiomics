.. _radiomics-removed-features-label:

==========================
Excluded Radiomic Features
==========================

Some commonly know features are not supported (anymore) in PyRadiomics. These features are listed here, so as to provide
a complete overview, as well as argumentation for why these features are excluded from PyRadiomics

Excluded GLCM Features
----------------------

For included features and class definition, see :ref:`radiomics-glcm-label`.

.. _radiomics-excluded-sumvariance-label:

Sum Variance
############

.. math::
    \textit{sum variance} = \displaystyle\sum^{2N_g}_{k=2}{(k-SA)^2p_{x+y}(k)}

Sum Variance is a measure of heterogeneity that places higher weights on
neighboring intensity level pairs that deviate more from the mean.

This feature has been removed, as it is mathematically identical to Cluster Tendency (see
:py:func:`~radiomics.glcm.RadiomicsGLCM.getClusterTendencyFeatureValue()`).

The mathematical proof is as follows:

(1) As defined in GLCM,
    :math:`p_{x+y}(k) = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)},\text{ where }i+j=k, k \in \{2, 3, \dots, 2N_g\}`

(2) Starting with cluster tendency as defined in GLCM:

.. math::
    \textit{cluster tendency} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}
      {\big(i+j-\mu_x-\mu_y\big)^2p(i,j)}

    = \displaystyle\sum^{2N_g}_{k=2}{\Big[\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}
      {\big(i+j-\mu_x-\mu_y\big)^2p(i,j)} \text{, where }i+j=k\Big]}

    = \displaystyle\sum^{2N_g}_{k=2}{\Big[\displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}
      {\big(k-(\mu_x+\mu_y)\big)^2p(i,j)} \text{, where }i+j=k \Big]}

.. note::
    Because inside the sum :math:`\sum^{2N_g}_{k=2}`, :math:`k` is a constant, and so are :math:`\mu_x` and
    :math:`\mu_y`, :math:`\big(k-(\mu_x+\mu_y)\big)^2` is constant and can be taken outside the inner sum
    :math:`\sum^{N_g}_{i=1}\sum^{N_g}_{j=1}`.

.. math::
    = \displaystyle\sum^{2N_g}_{k=2}{\Big[\big(k-(\mu_x+\mu_y)\big)^2\displaystyle\sum^{N_g}_{i=1}
      \displaystyle\sum^{N_g}_{j=1}{p(i,j)} \text{, where }i+j=k \Big]}

(3) Using (1.) and (2.)

.. math::
    \textit{cluster tendency} = \displaystyle\sum^{2N_g}_{k=2}{\Big[\big(k-(\mu_x+\mu_y)\big)^2p_{x+y}(k)\Big]}

(4) As defined in GLCM, :math:`p_x(i) = \sum^{N_g}_{j=1}{P(i,j)}` and :math:`\mu_x = \sum^{N_g}_{i=1}{p_x(i)i}`,
    therefore :math:`\mu_x = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{P(i,j)i}`

(5) Similarly as in (4.), :math:`\mu_y = \sum^{N_g}_{j=1}\sum^{N_g}_{i=1}{P(i,j)j}`

(6) Using (4.) and (5.), :math:`\mu_x` and :math:`\mu_y` can then be combined as follows:

.. math::
    \mu_x + \mu_y = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{P(i,j)i} +
    \displaystyle\sum^{N_g}_{j=1}\displaystyle\sum^{N_g}_{i=1}{P(i,j)j}

    = \displaystyle\sum^{N_g}_{j=1}\displaystyle\sum^{N_g}_{i=1}{P(i,j)i + P(i, j)j}

    = \displaystyle\sum^{N_g}_{j=1}\displaystyle\sum^{N_g}_{i=1}{P(i,j)(i + j)}

    = \displaystyle\sum^{2N_g}_{k=2}{\Big[\displaystyle\sum^{N_g}_{j=1}
    \displaystyle\sum^{N_g}_{i=1}{P(i,j)(i + j)} \text{, where } k = i + j\Big]}

    = \displaystyle\sum^{2N_g}_{k=2}{p_{x+y}(k)k} = \textit{sum average (SA)}

(7) Combining (3.) and (6.) yields the following formula:

.. math::
    \text{Cluster Tendency} =
    \displaystyle\sum^{2N_g}_{k=2}{\Big[\big(k-SA\big)^2p_{x+y}(k)\Big]} =
    \textit{ sum variance}

Q.E.D

.. _radiomics-excluded-dissimilarity-label:

Dissimilarity
#############

.. math::
    \textit{dissimilarity} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)|i-j|}

Dissimilarity is a measure of local intensity variation defined as the mean absolute difference between the
neighbouring pairs. A larger value correlates with a greater disparity in intensity values
among neighboring voxels.

This feature has been removed, as it is mathematically identical to Difference Average (see
:py:func:`~radiomics.glcm.RadiomicsGLCM.getDifferenceAverageFeatureValue()`).

The mathematical proof is as follows:

(1) As defined in GLCM, :math:`p_{x-y}(k) = \sum^{N_g}_{i=1}\sum^{N_g}_{j=1}{p(i,j)},\text{ where }|i-j|=k`

(2) Starting with Dissimilarity as defined in GLCM:

.. math::
    \textit{dissimilarity} = \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)|i-j|}

    = \displaystyle\sum^{N_g-1}_{k=0}{\Big[
    \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)|i-j|} \text{, where }|i-j|=k\Big]}

    = \displaystyle\sum^{N_g-1}_{k=0}{\Big[
    \displaystyle\sum^{N_g}_{i=1}\displaystyle\sum^{N_g}_{j=1}{p(i,j)k} \text{, where }|i-j|=k\Big]}

(3) Using (1.) and (2.)

.. math::
    \textit{dissimilarity} = \displaystyle\sum^{N_g-1}_{k=0}{p_{x-y}(k)k} = \textit{difference average}

Q.E.D.

.. _radiomics-excluded-gldm-label:

Excluded GLDM Features
----------------------

For included features and class definition, see :ref:`radiomics-gldm-label`.

.. _radiomics-excluded-gldm-dependence-percentage-label:

Dependence percentage
#####################

.. math::
    \textit{dependence percentage} = \frac{N_z}{N_p}

Dependence percentage is the ratio between voxels with a dependence zone and the total number of voxels in the image.
Because PyRadiomics allows for incomplete dependence zones, all voxels have a dependence zone and :math:`N_z = N_p`.
Therefore, this feature would always compute to 1.

.. _radiomics-excluded-gldm-glnn-label:

Gray Level Non-Uniformity Normalized
####################################

.. math::
    \textit{GLNN} = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_d}_{j=1}{\textbf{P}(i,j)}\right)^2}{N_z^2}

Measures the similarity of gray-level intensity values in the image, where a lower GLNN value
correlates with a greater similarity in intensity values. This is the normalized version of the GLN formula.

This formula has been removed, because due to the definition of GLDM matrix (allowing incomplete zones), this feature is
equal to first order Uniformity (see :py:func:`~radiomics.firstorder.RadiomicsFirstOrder.getUniformityFeatureValue()`).

The mathematical proof is as follows:

(1) Starting with Gray Level Non-Uniformity Normalized as defined in GLDM,

.. math::
    \textit{GLNN} = \frac{\sum^{N_g}_{i=1}\left(\sum^{N_d}_{j=1}{\textbf{P}(i,j)}\right)^2}{N_z^2}

    = \displaystyle\sum^{N_g}_{i=1}{
        \frac{ \left( \sum^{N_d}_{j=1}{ \textbf{P}(i,j) } \right)^2 }{ N_z^2 }
    }

    = \displaystyle\sum^{N_g}_{i=1}{ \left(
        \frac{ \sum^{N_d}_{j=1}{ \textbf{P}(i,j) } }{ N_z }
    \right)^2}

    = \displaystyle\sum^{N_g}_{i=1}{ \left(
        \sum^{N_d}_{j=1}{ \frac{ \textbf{P}(i,j) } { N_z } }
    \right)^2}

(2) As defined in GLDM, :math:`p(i,j) = \frac{\textbf{P}(i,j)}{N_z}`

(3) Using (1.) and (2.)

.. math::
    \textit{GLNN} = \displaystyle\sum^{N_g}_{i=1}{ \left(
        \sum^{N_d}_{j=1}{ p(i,j) }
    \right)^2}

(4) Because in the PyRadiomics definition incomplete dependence zones are allowed, every voxel in the ROI has a
    dependence zone. Therefore, :math:`N_z = N_p` and :math:`\sum^{N_d}_{j=1}{\textbf{P}(i,j)}` equals the number of voxels
    with gray level :math:`i` and is equal to :math:`\textbf{P}(i)`, the first order histogram with :math:`N_g` discrete
    gray levels, as defined in first order.

(5) As defined in first order, :math:`p(i) = \frac{\textbf{P}(i)}{N_p}`

(6) Using (2.), (4.) and (5.)

.. math::
    \displaystyle\sum^{N_d}_{j=1}{\textbf{P}(i,j)} = \textbf{P}(i)

    \frac{\sum^{N_d}_{j=1}{\textbf{P}(i,j)}}{N_z} = \frac{\textbf{P}(i)}{N_p}

    \displaystyle\sum^{N_d}_{j=1}{\frac{\textbf{P}(i,j)}{N_z}} = \frac{\textbf{P}(i)}{N_p}

    \displaystyle\sum^{N_d}_{j=1}{p(i,j)} = p(i)

(7) Combining (3.) and (6.) yields:

.. math::
    \textit{GLNN} = \displaystyle\sum^{N_g}_{i=1}{p(i)^2} = Uniformity

Q.E.D.
