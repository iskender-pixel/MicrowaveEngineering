import skrf as rf
import numpy as np

#------Variables-------
epseff = 1.87
l = 20e-3

#------Functions------

def Sparam_to_ABCD(S, Z0):


    # entries = len(S[:,:,0])
    ABCD = rf.network.s2a(S, Z0)


    return ABCD




def ReadTXT(Path):
    frequency = []
    real_part = []
    imaginary_part = []

    # Open the file and read line by line
    with open(Path, 'r') as file:
        # Skip the header lines (first two lines)
        next(file)  # Skip the first line (header)
        next(file)  # Skip the second line (separator)

        # Process each line in the file
        for line in file:
            # Skip lines that start with '#' or are empty
            if line.strip().startswith('#') or not line.strip():
                continue

            # Split the line into columns based on whitespace
            columns = line.split()

            # Ensure the line has exactly 3 columns (frequency, real, imaginary)
            if len(columns) == 2 or len(columns) == 3:
                # Append the values to the respective lists
                frequency.append(float(columns[0]))
                real_part.append(float(columns[1]))
                if len(columns) == 3:
                    imaginary_part.append(float(columns[2]))


    return np.array(frequency), np.array(real_part), np.array(imaginary_part)

def prepareSMatrix(s11Real, s11Imag, s12Real, s12Imag):

    entries = len(s11Real)
    S = np.zeros((entries,2,2),dtype=complex)

    s11 = s11Real + 1j * s11Imag
    s12 = s12Real + 1j * s12Imag
    s21 = s12
    s22 = s11

    for i in range(entries):
        Si = np.array([[s11[i], s12[i]],[s21[i], s22[i]]])
        S[i,:,:] = Si

    return S
if __name__ == '__main__':

    #Get the values for Z0
    _, Z0Real, _= ReadTXT("A1_Z0.txt")
    Z0 = Z0Real

    #Get the values for the S11
    f, S11Real, S11Imag = ReadTXT("A1_S11_parameters_RealImag.txt")
    _ , S12Real, S12Imag = ReadTXT("A1_S12_parameters_RealImag.txt")

    S = prepareSMatrix(S11Real, S11Imag, S12Real, S12Imag)

    ABCD = Sparam_to_ABCD(S, Z0)
    print(ABCD)