# Chito-oligosaccharide mass spectrometry - cosms

## Installation (conda)

Run following steps in Anaconda Prompt

1. Create new Python 3 environment

    ```
    conda create -n cosmspym2 python=3
    ```

2. Activate new environment

    ```
    conda activate cosmspym2
    ```

3. Install cosms and dependencies

    Go to directory containing *setup.py* and run

    ```
    pip install -e .
    ```

4. Install dependencies for Jupyter

    ```
    conda install -c anaconda ipykernel
    conda install -c anaconda jupyter
    ```

5. Register environment

    ```
    python -m ipykernel install --user --name=cosmspym2 --display-name "Python 3 (cosms pymzml2)"
    ```

## Run Jupyter with cosms

Run following steps in Anaconda Prompt

1. Activate environment

    ```
    conda activate cosmspym2
    ```

2. Run Jupyter notebook

    ```
    jupyter notebook
    ```