# import re
# import numpy as np
# import matplotlib.pyplot as plt

# def extract_data_from_file(filename):
#     # pattern = re.compile(r"total gas generation rate \[mol \* m-3\]:\s*(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
#     pattern0 = re.compile(r">>> Sac mass (* alpha): \s*(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
#     pattern1 = re.compile(r"Sac mass (w/o alpha): \s*(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
#     data0 = []
#     data1 = []

#     with open(filename, 'r') as file:
#         for line in file:
#             match0 = pattern0.search(line)
#             match1 = pattern1.search(line)
#             if match0:
#                 value = float(match0.group(1))
#                 data0.append(value)
#             if match1:
#                 value = float(match1.group(1))
#                 data1.append(value)

#     return np.array(data0), np.array(data1)

# def plot_data(data0, data1):
#     x = np.arange(len(data0))
#     plt.figure(figsize=(10, 6))
#     plt.plot(x, data0, marker='o', linestyle='-', color='b')
#     plt.plot(x, data1, marker='o', linestyle='-', color='r')
#     plt.title("Total Gas Generation Rate Over Time")
#     plt.xlabel("Sample Index")
#     plt.ylabel("Rate [mol * m^-3]")
#     plt.grid(True)
#     plt.tight_layout()
#     plt.show()

# if __name__ == "__main__":
#     filename = "log.interADMFoam"  # replace with your actual file name
#     data0, data1 = extract_data_from_file(filename)
#     if data.size > 0:
#         plot_data(data0[1::10], data1[1::10])
#     else:
#         print("No valid data found in the file.")


import re
import numpy as np
import matplotlib.pyplot as plt

def extract_data_from_file(filename):
    # Updated regex to handle literal asterisk and optional whitespace
    pattern0 = re.compile(r">>> Sac mass \(\* alpha\):\s*(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
    pattern1 = re.compile(r"Sac mass \(w/o alpha\):\s*(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")
    
    data0 = []
    data1 = []

    with open(filename, 'r') as file:
        for line in file:
            match0 = pattern0.search(line)
            match1 = pattern1.search(line)
            if match0:
                value = float(match0.group(1))
                data0.append(value)
            if match1:
                value = float(match1.group(1))
                data1.append(value)

    return np.array(data0), np.array(data1)

def plot_data(data0, data1):
    x = np.arange(len(data0))
    plt.figure(figsize=(10, 6))
    plt.plot(x, data0, marker='', linestyle='-', color='b', label='Sac mass (* alpha)')
    plt.plot(x, data1, marker='', linestyle='-', color='r', label='Sac mass (w/o alpha)')
    # plt.ylim(8.56e-8, 8.59e-8)
    plt.ylim(8.56e-5, 8.59e-5)
    plt.title("Sac Mass Comparison Over Time")
    plt.xlabel("Sample Index")
    plt.ylabel("Sac Mass")
    plt.legend(['Sac Mass with real time phase volume','Sac Mass with constant phase volume'])
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    filename = "log.interADMFoam"  # replace with your actual file name
    data0, data1 = extract_data_from_file(filename)

    if data0.size > 0 and data1.size > 0:
        # plot_data(data0[1::10], data1[1::10])  # Sample every 10th point starting from index 1
        plot_data(1e3*data0[1::10], 1e3*data1[1::10])  # Sample every 10th point starting from index 1
    else:
        print("No valid data found in the file.")
