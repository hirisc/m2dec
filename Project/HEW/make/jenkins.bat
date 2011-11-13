mkdir ..\libm2dec\Release
mkdir ..\libm2dec\Debug
mkdir ..\h264dec\Release
mkdir ..\h264dec\Debug
mkdir ..\m2dec\Release
mkdir ..\m2dec\Debug
"C:\Program Files (x86)\Renesas\Hew\hmake.exe" m2dec.mak > buildlog.txt
find "Make process completed" buildlog.txt
if errorlevel 1 exit /b 1
exit /b 0
