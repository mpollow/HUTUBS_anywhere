# HUTUBS_anywhere
This project allows to derive HRTF datasets from the HUTUBS database for any desired set of directions. It uses the HRIR representation given as SH-coefficients and converts them to a NumPy array for a given set of angles (theta, phi).

The HUTUBS HRIR datasets should be placed in the folder HUTUBS_HRIRs. Download the dataset HRIRs.zip from the [HUTUBS project website](https://depositonce.tu-berlin.de/items/dc2a3076-a291-417e-97f0-7697e332c960)

## Example usage:
Use the function get_HRTFs_fromSH to get frequency-domain data for the requested directions. Here an example for getting HRIRs for the front direction, 45 degree to the left side and 90 degree to the left side.

```
import numpy as np
from HUTUBS_anywhere import get_HRTFs_fromSH

# angles given in rad
theta = np.array([np.pi/2, np.pi/2, np.pi/2])
phi = np.array([0, np.pi/4, np.pi/2])
HRTFs, sampling_rate = get_HRTFs_fromSH('./HUTUBS_HRIRs_2/pp21_SHcoefficients_simulated.mat', theta, phi)
HRIRs = np.fft.irfft(HRTFs)
```

