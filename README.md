<h3 align="center">ProtoSim</h3>

---

<p align="center"> - Calculates and visualizes the variation of internal protocell material over time.
- Implements specific differential equations to simulate protocell behavior.
- Allows users to customize input parameters and conditions.
    <br> 
</p>

## 🏁 Getting Started <a name = "getting_started"></a>

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### 💻 Download Repo

If you don't have Git installed on your system, follow these steps to install it:

#### Windows 

1. Download Git from the official website: [https://git-scm.com/](https://git-scm.com/).

2. Run the downloaded file (.exe) and follow the installer instructions.

3. During installation, make sure to select "Use Git from the Windows Command Prompt" to add Git to your PATH.

4. Verify the installation by running the following command in the command prompt or PowerShell:

    ```bash
    git --version
    ```
#### MacOS

1. Open the terminal.

2. Install Git using Homebrew (if you don't have Homebrew, follow the instructions at [https://brew.sh/](https://brew.sh/)):

    ```bash
    brew install git
    ```

3. Verify the installation by running the following command in the terminal:

    ```bash
    git --version
    ```
#### Linux (Debian/Ubuntu):

1. Open the terminal.

2. Run the following commands:

    ```bash
    sudo apt-get update
    sudo apt-get install git
    ```

3. Verify the installation by running the following command in the terminal:

    ```bash
    git --version
    ```

Now you are ready! Run the following command to clone repo:

```
git clone https://github.com/Piltxi/ProtoSim
```

### 🔮 Start simulation

```
cd ProtoSim/src
python3 main.py
```

### Note

📤 To enable printing of some additional information or debug printing: 

```
cd src
python3 main.py -v
```

📥 To delete the directory with output information: 

```
cd src
python3 main.py -r
```

ProtoSim reads the input parameters and reactions from two text files. 
The files are located in input/ .

## ✍️ Authors <a name = "authors"></a>

- [@piltxi](https://github.com/Piltxi/)
