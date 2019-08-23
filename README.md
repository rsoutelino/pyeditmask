pyeditmask
==========

GUI to edit land mask in ROMS grid

Quick start
----------
- install dependencies
  - wxpython, numpy, scipy, matplotlib, netCDF4

- or, clone this repo into /source and use the docker image:
```
docker run -ti --rm -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $HOME:$HOME -v /home/roms:/home/roms  -v /static:/static rsoutelino/pyeditmask
```

- Note that is necessary to mount `$HOME` because that's how the X forwarding is handled from host to container. 

- Note that a /static (or any other name of your preference) folder is also needed, so the grid file can be edited and persist once the docker container exits. 

- run `python pyeditmask.py`
- load the grid
- Using the "A" buttons: toggle the buttom and start clicking with the left mouse buttom to make roughly a closed poligon. Click the right buttom when you're finished.   
