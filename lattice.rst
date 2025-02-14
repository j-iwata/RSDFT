.. _lattice:

==============================
格子ベクトルと逆格子ベクトル
==============================

------------------------------------
定義と性質
------------------------------------

結晶の周期性および単位胞（Unit cell）を定める３つのベクトル

.. math::
  \boldsymbol{a}_1 =
  \begin{pmatrix}
  a_{1}^x \\
  a_{1}^y \\
  a_{1}^z
  \end{pmatrix}
  , \quad
  \boldsymbol{a}_2 =
  \begin{pmatrix}
  a_{2}^x \\
  a_{2}^y \\
  a_{2}^z
  \end{pmatrix}
  , \quad
  \boldsymbol{a}_3 =
  \begin{pmatrix}
  a_{3}^x \\
  a_{3}^y \\
  a_{3}^z
  \end{pmatrix}
  :label: avec

これを格子ベクトルと呼ぶ。格子ベクトルが実空間での周期性を表すのに対し、次で定義される逆格子ベクトル

.. math::
  \boldsymbol{b}_1 = 2\pi\frac{\boldsymbol{a}_2\times\boldsymbol{a}_3}
  {\boldsymbol{a}_1\cdot \left(\boldsymbol{a}_2\times\boldsymbol{a}_3\right)}
  , \quad
  \boldsymbol{b}_2 = 2\pi\frac{\boldsymbol{a}_3\times\boldsymbol{a}_1}
  {\boldsymbol{a}_1\cdot \left(\boldsymbol{a}_2\times\boldsymbol{a}_3\right)}
  , \quad
  \boldsymbol{b}_3 = 2\pi\frac{\boldsymbol{a}_1\times\boldsymbol{a}_2}
  {\boldsymbol{a}_1\cdot \left(\boldsymbol{a}_2\times\boldsymbol{a}_3\right)}
  :label: bvec

は、（フーリエ変換で移った先の）波数空間での周期性を表す。

格子ベクトルと逆格子ベクトルは、次の直交関係を満たす。

.. math::
  \boldsymbol{a}_i \cdot \boldsymbol{b}_j = 2\pi\delta_{ij}
  :label: ortho

いま、格子ベクトルと逆格子ベクトルの成分をそれぞれ並べ、以下のような行列を作ってみる
（これはカルテシアン座標と\ :ref:`格子座標<coordinates>`\ の間の座標変換の行列となっている）。

.. math::
  A \equiv
  \begin{pmatrix}
  a_{1}^x & a_{2}^x & a_{3}^x \\
  a_{1}^y & a_{2}^y & a_{3}^y \\
  a_{1}^z & a_{2}^z & a_{3}^z \\
  \end{pmatrix}
  , \quad
  B \equiv
  \begin{pmatrix}
  b_{1}^x & b_{2}^x & b_{3}^x \\
  b_{1}^y & b_{2}^y & b_{3}^y \\
  b_{1}^z & b_{2}^z & b_{3}^z \\
  \end{pmatrix}
  :label: abmat

式 :eq:`ortho` の直交性から、行列の積として

.. math::
  A^\top B
  =
  \begin{pmatrix}
  a_{1}^x & a_{1}^y & a_{1}^z \\
  a_{2}^x & a_{2}^y & a_{2}^z \\
  a_{3}^x & a_{3}^y & a_{3}^z \\
  \end{pmatrix}
  \begin{pmatrix}
  b_{1}^x & b_{2}^x & b_{3}^x \\
  b_{1}^y & b_{2}^y & b_{3}^y \\
  b_{1}^z & b_{2}^z & b_{3}^z \\
  \end{pmatrix}
  = 2\pi I
  :label: abortho

となることがわかる。これから、\ :math:`A, B` の逆行列は、互いの転置行列に一致することがわかる。

.. math::
  A^{-1} = \frac{1}{2\pi}B^\top
  , \quad
  B^{-1} = \frac{1}{2\pi}A^\top
  :label: inverse


------------------------------------
格子ベクトルの独立な成分
------------------------------------

式 :eq:`avec` は9つの成分を持っている。これは、何らかの計算を行ったりするときに必要となるが、
たんに物理系を記述するという目的だけであれば、系全体の並進と回転の自由度を取り除いた6つの成分で十分である。
例えば結晶構造データの国際的なフォーマットであるCIF（Crystallographic Information File）では

.. code-block::

  _cell_length_a   3.84928
  _cell_length_b   3.84928
  _cell_length_c   3.84928
  _cell_angle_alpha   60.00000
  _cell_angle_beta    60.00000
  _cell_angle_gamma   60.00000

のように、格子ベクトルの長さと、それぞれのベクトルのなす角度、という情報だけが記録されている。式 :eq:`avec` との対応を明示すれば

.. math::
  a = \|\boldsymbol{a}_1\|
  , \quad
  b = \|\boldsymbol{a}_2\|
  , \quad
  c = \|\boldsymbol{a}_3\|

.. math::
  \alpha = \arccos\left(\frac{\boldsymbol{a}_2\cdot\boldsymbol{a}_3}{\|\boldsymbol{a}_2\|\|\boldsymbol{a}_3\|}\right)
  , \quad
  \beta = \arccos\left(\frac{\boldsymbol{a}_3\cdot\boldsymbol{a}_1}{\|\boldsymbol{a}_3\|\|\boldsymbol{a}_1\|}\right)
  , \quad
  \gamma = \arccos\left(\frac{\boldsymbol{a}_1\cdot\boldsymbol{a}_2}{\|\boldsymbol{a}_1\|\|\boldsymbol{a}_2\|}\right)
  :label: len_and_ang
  
となる。
