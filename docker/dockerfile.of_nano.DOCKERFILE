FROM opencfd/openfoam-default:latest

# Install nano and clean up
RUN apt-get update && apt-get install -y --no-install-recommends nano \
    && rm -rf /var/lib/apt/lists/*

