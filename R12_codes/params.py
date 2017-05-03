import os
import socket

hostname = socket.gethostname()

EXT = '.png'

if hostname.endswith('astro.washington.edu'):
    TABLE_DIR = '/astro/users/philrose/python/tables/'
    PHAT_DIR_EXTRA = '/astro/net/angst2/philrose/PHAT/paper/supporting_files/'
    PHAT_CODE_DIR = '/astro/net/angst2/philrose/PHAT/code/'
    PHAT_TRACKS_DIR = PHAT_DIR_EXTRA
    PHAT_DIR_EXTRA2 = PHAT_DIR_EXTRA
    PHAT_TRACKS_BASE = '/astro/net/angst2/philrose/PHAT/tracks/'
    VW_BASE = PHAT_TRACKS_BASE
    TRACKS_BASE ='/astro/net/angst2/philrose/tracks/'
    PHAT_DATA_BASE = '/astro/net/angst2/philrose/PHAT/FullBricks/'
    PHAT_DIR_EXTRA3 = PHAT_DIR_EXTRA
else:
    PHAT_DIR_EXTRA = '/Users/phil/research/PHAT/'
    PHAT_CODE_DIR = '/Users/phil/research/PHAT/code/'
    PHAT_TRACKS_DIR = '/Users/phil/research/PHAT/UVBrightpaper/fromLeo/'
    PHAT_DIR_EXTRA2 = '/Users/phil/research/PHAT/UVBrightpaper/fromCliff/'
    PHAT_TRACKS_BASE = '/Users/phil/research/PHAT/tracks/'
    VW_BASE = '/Users/phil/research/PHAT/UVBrightpaper/fromLeo/'
    TRACKS_BASE ='/Users/phil/research/Italy/tracks/'
    PHAT_DATA_BASE = '/Users/phil/research/PHAT/Data/'
    PHAT_DIR_EXTRA3 = '/Users/phil/research/PHAT/Bulge/'
