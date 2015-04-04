if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:$PWD/lib
else
  setenv LD_LIBRARY_PATH $PWD/lib
endif
