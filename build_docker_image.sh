docker build -t mmpbsa .
docker run -v `pwd`:/workspace:Z -it mmpbsa bash
