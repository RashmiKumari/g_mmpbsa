docker build -t mmpbsa .
docker run -v ../:/workspace:Z -it mmpbsa bash
