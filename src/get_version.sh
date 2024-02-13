gID=`git log -n 1 | grep commit | cut -f 2 -d' ' | tr -d $'\n'` 
echo "      character(len=40),parameter,public :: Ash3d_GitComID ='${gID}'" > Ash3d_version.h 
