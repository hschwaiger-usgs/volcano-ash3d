echo -e "      character(len=40),parameter,public :: Ash3d_GitComID ='\c" > Ash3d_version.h
git log -n 1 | grep commit | cut -f 2 -d' ' | tr -d $'\n' >> Ash3d_version.h
echo -e "'\c" >> Ash3d_version.h
