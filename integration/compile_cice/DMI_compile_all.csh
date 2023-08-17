#!/bin/csh -f

foreach compiler (cray intel gnu)
foreach debug (true false) 
./Make_DMI_CICE.csh ${compiler} ${debug}
end
end
