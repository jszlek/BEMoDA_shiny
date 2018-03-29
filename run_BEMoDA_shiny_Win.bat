rem In order to use this *.bat file you need to set PATH to recognize the R command as a command-line tool.
rem 
rem  You can follow one of the methods provided by superuser.com community
rem  (https://superuser.com/questions/949560/how-do-i-set-system-environment-variables-in-windows-10)
rem  or try to
rem  click start > right click on "My Computer" > Choose "Properties" > "Advanced" >
rem  "Envonment variables" > Under system variables choose "Path" and click
rem  Edit > on variable value at the end of it add ;C:\Program Files\R\R-3.4.3\bin
rem  
rem  Click OK. And you are done.
rem  
rem  The C:\Program Files\R\R-3.4.3\bin is the place you have the R
rem  executables which you can find from explorer.
rem  
rem  After that your can call R from command prompt like
rem  Rscript.exe myRscript.R e.g


Rscript.exe run_BEMoDA_shiny.R

exit 0
