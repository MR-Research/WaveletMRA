# Wavelet-MRA

Discrete Wavelet Transform MRA analysis and synthesis written in MATLAB


## Installation

 - Clone the repository.
```bash
git clone https://github.com/MrRobot-oss/WaveletMRA
```
 - Open using MATLAB

## Usage

```matlab
dwt_mra(base_a, base_d, signal, level, maxlevel, support, res)
idwt_mra_v2(base_a, base_d, signal, level, maxlevel, support)
```
## Parameters

 1. base_a: approximation coefficients used for decomposition/reconstruction.
 2. base_d: detail coefficients used for decomposition/reconstruction.
 3. signal: 1-D signal to be analyzed.
 4. level: starting decomposition/analysis level, should always be 1
 5. maxlevel: last decomposition/analysis level, can be 1 up to base2log(N), where N is the length of signal.
 6. support: number of approximation and detail coefficients (only orthogonal wavelets are supported)
 7. res: matrix of size maxlevel x length(signal)/2. Used to store the decomposition approximations for each analyzed level in the recursive procedure.
 

## Limitations
The following are the current limitations of the program:
 - Only orthogonal wavelets supported.
 - Signal extension assumes a periodic signal.

## Known issues
No known issues have been reported so far. If you encounter an issue, please post a new issue so it can be verified and fixed.

## TODO

 - Add version that supports coefficient filter passed as argument

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## License
This program is distributed under the GNU-GPL version 3, you can read more about it here: [GNU-GPL](https://www.gnu.org/licenses/gpl-3.0.txt).