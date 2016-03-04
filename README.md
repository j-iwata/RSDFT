# RSDFT ( ver.1.2.1 )

Real-space finite-difference pseudopotential density-functional-theory code for first-principles material simulations on massively-parallel computers

# DIRECTORIES
```
doc/		User's guide and Tutorial (日本語版あり)  
examples/	Examples of RSDFT for solids (RSSOL) and molecules (RSMOL)  
src/		Source programs of RSDFT  
utility/	Utility programs  
```
# LICENSE

   Copyright 2015 Junichi Iwata

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

# REFERENCES

"A massively-parallel electronic-structure calculations based on real-space density functional theory"  
Jun-Ichi Iwata, Daisuke Takahashi, Atsushi Oshiyama, Taisuke Boku, Kenji Shiraishi, Susumu Okada, Kazuhiro Yabana  
Journal of Computational Physics 229, pp2339-2363 (2010).  


"Performace evaluation of ultra-large-scale first-principles electronic structure calculation code on the K computer"  
Yukihiro Hasegawa, Jun-Ichi Iwata, Miwako Tsuji, Daisuke Takahashi, Atsushi Oshiyama, Kazuo Minami, Taisuke Boku, Hikaru Inoue, Yoshito Kitazawa, Ikuo Miyoshi, Mitsuo Yokokawa  
Internationa Journal of High Performance Computing Applications 28, pp335-355, (2014).  


# ACKNOWLEDGEMENT

## Contributors
```
内田和之        Kazuyuki Uchida  
小泉健一        Ken-ichi Koizumi  
澤田啓介        Keisuke Sawada  
重田育照        Yasuteru Shigeta  
高橋大介        Daisuke Takahashi  
辻美和子        Miwako Tsuji  
長谷川幸弘      Yukihiro Hasegawa  
古家真之介      Shinnosuke Furuya  
李 瀚           Han Li  
矢花一浩        Kazuhiro Yabana  
押山 淳         Atsushi Oshiyama  
```
## MINPACK ( src/ext1/ )

This product includes software developed by the University of Chicago, as Operator of Argonne National Laboratory.

## FFTE ( src/ext2/ )

This product includes a fast Fourier transform package FFTE
( http://www.ffte.jp )


