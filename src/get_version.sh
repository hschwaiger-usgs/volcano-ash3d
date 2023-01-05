echo -n "      character(len=40),parameter,public :: Ash3d_GitComID ='" > Ash3d_version.h
git log -n 1 | grep commit | cut -f 2 -d' ' | tr -d $'\n' >> Ash3d_version.h
echo -n "'" >> Ash3d_version.h
