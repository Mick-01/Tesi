FROM opencfd/openfoam-run:latest

# Install nano and ParaView
RUN apt-get update && apt-get install -y --no-install-recommends \
      nano \
      paraview \
    && rm -rf /var/lib/apt/lists/*

