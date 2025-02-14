.. _coordinates:

==============================
格子座標
==============================

------------------------------------
XYZ座標と格子座標
------------------------------------

格子座標は、\ :ref:`格子ベクトル<lattice>`\ を基底にとった位置の指定の仕方である。カルテシアン座標（デカルト座標またはXYZ座標）で

.. math::
  \boldsymbol{r}
  =
  \begin{pmatrix}
  x \\
  y \\
  z
  \end{pmatrix}
  :label: rpos

と表される位置は、格子座標では

.. math::
  \boldsymbol{r}
  =
  \alpha_1 \boldsymbol{a}_1 +
  \alpha_2 \boldsymbol{a}_2 +
  \alpha_3 \boldsymbol{a}_3
  :label: rlat

と表される。 :math:`\boldsymbol{r}` が単位胞内の点であれば、\ :math:`\alpha_i` は

.. math::
  0 \leq \alpha_i \leq 1
  :label: alpha

という範囲の値をとることになる
（ただし :math:`\alpha_i = 0` と :math:`\alpha_i = 1` は周期境界条件により等価な点を指すことに注意する）。

.. .. math::
..   \begin{pmatrix}
..   x_1 \\
..   x_2 \\
..   x_3
..   \end{pmatrix}
..   =
..   \begin{pmatrix}
..   a_{11} & a_{12} & a_{13} \\
..   a_{21} & a_{22} & a_{23} \\
..   a_{31} & a_{32} & a_{33} \\
..   \end{pmatrix}
..   \begin{pmatrix}
..   \alpha_1  \\
..   \alpha_2 \\
..   \alpha_3
..   \end{pmatrix}

.. 座標変換の行列 :math:`A` は、格子ベクトルのxyz成分

.. .. math::
..   \boldsymbol{a}_1 =
..   \begin{pmatrix}
..   a_{1}^x \\
..   a_{1}^y \\
..   a_{1}^z
..   \end{pmatrix}
..   \equiv
..   \begin{pmatrix}
..   a_{11} \\
..   a_{21} \\
..   a_{31}
..   \end{pmatrix}
..   , \quad
..   \boldsymbol{a}_2 =
..   \begin{pmatrix}
..   a_{2}^x \\
..   a_{2}^y \\
..   a_{2}^z
..   \end{pmatrix}
..   \equiv
..   \begin{pmatrix}
..   a_{12} \\
..   a_{22} \\
..   a_{32}
..   \end{pmatrix}
..   , \quad
..   \boldsymbol{a}_3 =
..   \begin{pmatrix}
..   a_{3}^x \\
..   a_{3}^y \\
..   a_{3}^z
..   \end{pmatrix}
..   \equiv
..   \begin{pmatrix}
..   a_{13} \\
..   a_{23} \\
..   a_{33}
..   \end{pmatrix}

.. を並べたものである。同様に逆格子ベクトルのxyz成分を並べた行列を考え、

.. .. math::
..   B \equiv
..   \begin{pmatrix}
..   b_{11} & b_{12} & b_{13} \\
..   b_{21} & b_{22} & b_{23} \\
..   b_{31} & b_{32} & b_{33} \\
..   \end{pmatrix}

.. とおく。格子ベクトルと逆格子ベクトルは

.. .. math::
..   \boldsymbol{a}_i \cdot \boldsymbol{b}_j = 2\pi\delta_{ij}

.. という関係を満たすので

.. .. math::
..   \begin{pmatrix}
..   a_{11} & a_{21} & a_{31} \\
..   a_{12} & a_{22} & a_{32} \\
..   a_{13} & a_{23} & a_{33} \\
..   \end{pmatrix}
..   \begin{pmatrix}
..   b_{11} & b_{12} & b_{13} \\
..   b_{21} & b_{22} & b_{23} \\
..   b_{31} & b_{32} & b_{33} \\
..   \end{pmatrix}
..   =
..   A^\top B = 2\pi I

.. .. math::
..   \frac{1}{2\pi}A^\top = B^{-1}
..   \\
..   \frac{1}{2\pi}B^\top = A^{-1}

.. のように、転置行列が互いの逆行列となることがわかる。

.. RSDFTのプログラム内部では、原子座標やグリッド点は基本的に格子座標で扱っている。
.. ただし基底を :math:`\boldsymbol{a}_1, \boldsymbol{a}_2, \boldsymbol{a}_3` 
.. そのものではなく、各方向の単位ベクトルにとっている場合もある。

.. ------------------------------------
.. 微分演算子の変換
.. ------------------------------------

.. RSDFTでは微分演算子に対して次の座標変換

.. .. math::
..   \frac{\partial f}{\partial x_i}
..   = \sum_{j=1}^3 \frac{\partial f}{\partial \alpha_j}\frac{\partial \alpha_j}{\partial x_i}
..   = \sum_{j=1}^3 \frac{\partial f}{\partial \alpha_j} \frac{b_{ij}}{2\pi}

.. を行なった後 :math:`\frac{\partial f}{\partial \alpha_i}` 等に対して差分化を行う。２階微分は

.. .. math::
..   \frac{\partial^2 f}{\partial x_i^2}
..   &=
..   \frac{\partial}{\partial x_i}
..   \left(
..   \sum_{j=1}^3 \frac{\partial f}{\partial \alpha_j} \frac{b_{ij}}{2\pi}
..   \right) \\
..   &=
..   \sum_{k=1}^3 \frac{\partial}{\partial \alpha_k}
..   \left(
..   \sum_{j=1}^3 \frac{\partial f}{\partial \alpha_j} \frac{b_{ij}}{2\pi}
..   \right)
..   \frac{b_{ik}}{2\pi} \\
..   &=
..   \sum_{k=1}^3
..   \sum_{j=1}^3 \frac{\partial^2 f}{\partial \alpha_j \partial \alpha_k} \frac{b_{ij}}{2\pi}
..   \frac{b_{ik}}{2\pi} \\

.. となり、格子座標が斜交座標の場合クロスターム
.. （ :math:`\frac{\partial^2 f}{\partial \alpha_1 \partial \alpha_2}`
.. 等 ）が現れる。
