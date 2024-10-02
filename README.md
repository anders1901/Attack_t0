# Uncompressing Dilithium's public key

This repository contains the artifact related to the corresponding article.

## Prerequisites

The following material requires:

* Python >= 3.0
* Pipenv
* lpsolve5.5
* Dilithium (commit 3e9b9f1)

### Installation Steps
1. Check Python version:  
```bash
python --version
```
If not installed or outdated, you can install it with:
```bash
sudo apt-get install python3-full
```

2. Check Pip version:  
```bash
pip --version
```
If not installed, use:
```bash
sudo apt install python3-pip
```

3. Check Pipenv version:  
```bash
pipenv --version
```
If not installed, run:
```bash
sudo pip install --user pipenv --break-system-packages  
```   
If a Warning appears saying that pipenv and pipenv-resolver are installed somewhere not in PATH you can run:
```bash
echo 'export PATH="<path_to_local_pipenv_bin>:$PATH"' >> ~/.bash_profile  
source ~/.bash_profile
```

4. Install lpsolve with the provided package (on linux):
```bash
cd lp_solve_5.5/lpsolve55
sh ccc
```
This will create a shared library file in ```lp_solve_5.5/lpsolve55/bin/<rest_of_path>```


5. The Dilithium reference implementation requires OpenSSL, install it with:
```bash
sudo apt-get install openssl libssl-dev
```


## Uncompressing t0 from standard signatures

1. Install Virtual Environment Dependencies in `Attack_t0/`:
```bash 
pipenv install
```

2. Activate the Virtual Environment: 
```bash 
pipenv shell
```

3. Add Shared Library to Dynamic Library Path:
```bash 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/lp_solve_5.5/lpsolve55/bin/<rest_of_path>
```

4. Launch Jupyter Notebook 
```bash 
jupyter notebook
```

5. Go to ```Additional_functions/notebook/```,  open the file ```Attack_t0.ipynb``` and follow the instructions.


## Files

| Name                   | Description                                              |
| :---                   | :---                                                     |
| `Additional_functions` | Additional material used for the attack                  |
| `dilithium-master`     | Reference implementation of Dilithium from   [GitHub](https://github.com/pq-crystals/dilithium)               |
| `lp_solve_5.5`         | The library used to solve LP systems                     |


## License

This work is licensed under a [Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).

[![CC BY 4.0](https://i.creativecommons.org/l/by/4.0/88x31.png)](http://creativecommons.org/licenses/by/4.0/)

See [LICENSE.txt](./LICENSE.txt).

This artifact uses the Dilithium reference implementation from [GitHub](https://github.com/pq-crystals/dilithium), under the Apache 2.0 License, as a submodule.

