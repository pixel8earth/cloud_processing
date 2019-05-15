docker run --rm -v `pwd`/scripts/:/scripts -v `pwd`/streams:/streams pixel8earth/insanecloudposse:latest python /scripts/process_clouds.py $1 
