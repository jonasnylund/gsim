from collections.abc import Sequence
import time

import numpy as np
import pygame

from gsim.model import py_model  # type: ignore


class Viewer:
  """Viewer of particle clouds."""

  def __init__(self):
    self.display = pygame.display.set_mode(
      size=(1000, 500),
      flags=pygame.RESIZABLE,
    )
    self.rendering_time: float = 0

  def __del__(self):
    pygame.quit()

  def show(self, particles: Sequence[py_model.Particle]) -> None:
    """Renders all the given particles to the display."""
    start_time = time.time()
    size = self.display.get_size()
    center_x = size[0] / 2
    center_y = size[1] / 2

    self.display.fill((0, 0, 0))
    for particle in particles:
      x, y = particle.position[:2]
      pygame.draw.circle(
        self.display,
        color=(200, 200, 255),
        center=(x + center_x, y + center_y),
        radius=np.ceil(np.sqrt(particle.mass)),
      )
    pygame.display.flip()

    self.rendering_time += time.time() - start_time

  def update(self) -> bool:
    """Updates the display and polls the events."""
    stop = False
    for event in pygame.event.get():
      if event.type == pygame.QUIT:
        stop = True

    pygame.event.clear()
    return not stop