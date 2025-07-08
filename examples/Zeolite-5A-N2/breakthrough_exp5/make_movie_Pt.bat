del column_movie_Pt.mp4
set /A argVec[1]=1
set /A argVec[2]=1200
set /A argVec[3]=800
set /A argVec[4]=18
setlocal enabledelayedexpansion
set argCount=0
for %%x in (%*) do (
   set /A argCount+=1
   set "argVec[!argCount!]=%%~x"'n)
set PATH=%PATH%;C:\Program Files\gnuplot\bin;C:\Program Files\ffmpeg-master-latest-win64-gpl\bin;C:\Program Files\ffmpeg\bin
gnuplot.exe -c plot_column_Pt %argVec[1]% %argVec[2]% %argVec[3]% | ffmpeg.exe -f png_pipe -s:v "%argVec[2]%,%argVec[3]%" -i pipe: -c:v libx264 -pix_fmt yuv420p -crf %argVec[4]% -c:a aac column_movie_Pt.mp4
