docker build -t mmpbsa -f Dockerfile
docker run -v ../:/workspace:Z -it mmpbsa bash
