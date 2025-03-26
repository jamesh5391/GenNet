FROM continuumio/miniconda3

WORKDIR /workspace

COPY requirements_GenNet.txt .

# Create the environment
RUN conda create -n env_GenNet python=3.10.12 -y

# Install pip packages
RUN /bin/bash -c "source activate env_GenNet && \
    pip install --upgrade pip && \
    pip install -r requirements_GenNet.txt"

COPY . /workspace

# RUN mkdir -p /workspace/processed_data /workspace/results

RUN echo "conda activate env_GenNet" >> ~/.bashrc

# Set CMD to launch bash with environment activated
CMD ["/bin/bash", "-c", "source activate env_GenNet && exec bash"]
