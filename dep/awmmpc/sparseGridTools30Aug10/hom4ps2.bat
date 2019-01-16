
@bin\sym2num.exe < %1 > input.num

@move input.num bin\input.num

@bin\flwcrv.exe
