podman build -t mmpbsa -f Dockerfile
podman run -v ../:/workspace:Z -it mmpbsa bash
