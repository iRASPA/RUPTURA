cl /W4 /EHsc /Zi /std:c++20 /O2 /Zc:__cplusplus ^
/I "C:\msys64\mingw64\lib\python3.12\site-packages\pybind11\include" ^
/I "C:\msys64\mingw64\include\python3.12" ^
*.cpp ^
/link /out:ruptura.exe
pause

