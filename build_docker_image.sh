docker build -t mmpbsa -f Dockerfile
docker run -v `pwd`:/workspace:Z -it mmpbsa bash
