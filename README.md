DBA
===

Averaging for Dynamic Time Warping

This repository gives you different versions of DBA for different programming language, whether you want to have a warping window or not, etc. 

Each file corresponds to one of these combinations; if one is missing for your usage, let me know. Currently the length is limited to 1,000 (let me know if you need more). 

* `DBA.java` standard DBA in Java with no warping window and memory allocated statically
* `DBAWarpingWindow.java` same as `DBA.java` but with a warping window as a parameter
* `DBAMultiThreadCapable.java` same as `DBA.java` (still with no windows) but with memory allocated for each call to DBA. It can thus be called by several threads at the same time. 
* `DBAWarpingWindowMultiThreadCapable.java` same as `DBAWarpingWindow.java` (with warping window) but with memory allocated for each call to DBA. It can thus be called by several threads at the same time. 
* `DBA.m` Matlab implementation of DBA with no windows 
* `DBA.py` Python implementation of DBA with no windows
* `DBA_DTW.py` Python implementation of DBA with warping window
