import sys
from math import *
from subprocess import *
from pathlib import Path
from shutil import rmtree

n_frames = 1024
resolution = "720p"

resolutions = {
    "2160p": (3840, 2160),
    "1440p": (2560, 1440),
    "1080p": (1920, 1080),
    "720p": (1280, 720),
    "480p": (854, 480),
    "360p": (640, 360),
    "240p": (426, 240),
}

width, height = resolutions[resolution]

basepath = Path("frames")
rmtree(basepath, ignore_errors=True)
basepath.mkdir(exist_ok=True)

def render_frames(indices):
    ps = []
    for i in indices:
        print("Launching", i + 1, "/", n_frames)
        filepath = basepath / "frame_{:05d}.png".format(i)
        t = i * 2*pi / n_frames
        p1 = Popen([
            "./bin/render",
            "--time", str(t),
            "--width", str(width),
            "--height", str(height)],
            stdout=PIPE, stderr=DEVNULL)
        p2 = Popen(["convert", "pgm:-", str(filepath)], stdin=p1.stdout)
        ps.append(p2)

    for i, p in zip(indices, ps):
        print("Waiting for", i + 1, "/", len(ps))
        p.communicate()

n_processes = 14
n_rendered = 0
while n_rendered < n_frames:
    end = min(n_rendered + n_processes, n_frames)
    render_frames(list(range(n_rendered, end)))
    n_rendered += n_processes

p = Popen([
    "ffmpeg", "-framerate", "30",
    "-i", str(basepath / "frame_%05d.png"),
    "-s:v", "{}x{}".format(width, height),
    "-c:v", "libx264",
    "-profile:v", "high",
    "-crf", "20",
    "-pix_fmt", "yuv420p",
    "-y",
    "output.mp4"])
p.communicate()
