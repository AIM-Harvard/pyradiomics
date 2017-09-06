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

1. Sum Variance
###############

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
