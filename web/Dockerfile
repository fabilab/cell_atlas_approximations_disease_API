# Use the official Python image from the Docker Hub
FROM python:3.11-slim

# Set environment variables to avoid interactive prompts during package installation
ENV PYTHONUNBUFFERED=1

# By default, listen on port 5000
EXPOSE 5000/tcp

# Set the working directory in the container
WORKDIR /app

# Copy the requirements.txt file into the container at /app
COPY requirements.txt .
COPY secret_key.txt .

# Install the dependencies
RUN pip install -r requirements.txt

# Static assets
ADD static static

# Create and copy manifest
COPY create_manifest.py .
RUN python create_manifest.py

# Modules for backend
ADD api api
ADD models models

# App file
COPY config.yml .
COPY config.py .
COPY app.py .

# Specify the command to run on container start
CMD [ "python", "./app.py" ]
