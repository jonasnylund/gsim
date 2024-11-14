from collections.abc import Sequence
import time

import numpy as np
import pygame

from gsim.model import py_model  # type: ignore


def draw_circle_alpha(surface, color, center, radius):
  target_rect = pygame.Rect(center, (0, 0)).inflate((radius * 2 + 1, radius * 2 + 1))
  intermediate = pygame.Surface(target_rect.size, pygame.SRCALPHA)
  pygame.draw.circle(intermediate, color, (radius + 1, radius + 1), radius)
  surface.blit(intermediate, target_rect)


class Viewer:
  """Viewer of particle clouds."""

  def __init__(self):
    self.display = pygame.display.set_mode(
      size=(1000, 500),
      flags=pygame.RESIZABLE | pygame.DOUBLEBUF,
    )
    self.scale: float = 1
    self.center: tuple[float, float] = (0, 0)
    self.rendering_time: float = 0

  def __del__(self):
    pygame.quit()

  def show(self, particles: Sequence[py_model.Particle]) -> None:
    """Renders all the given particles to the display."""
    start_time = time.time()
    size = self.display.get_size()
    center_x = size[0] / 2 - self.center[0]
    center_y = size[1] / 2 - self.center[1]

    self.display.fill((0, 0, 0))
    for particle in particles:
      x, y = particle.position[:2]
      size = np.pow(particle.mass, 1 / 3) * 2 * self.scale
      opacity = min(size, 2) / 2 * 255

      draw_circle_alpha(
        self.display,
        color=(200, 200, 255, int(opacity)),
        center=(x * self.scale + center_x, y * self.scale + center_y),
        radius=max(1, size),
      )
    pygame.display.flip()

    self.rendering_time += time.time() - start_time

  def update(self) -> bool:
    """Updates the display and polls the events."""
    stop = False
    for event in pygame.event.get():
      if event.type == pygame.QUIT:
        stop = True
      elif event.type == pygame.MOUSEWHEEL:
        self.scale *= 1 + np.clip(event.y, -2, 2) / 10
      elif event.type == pygame.MOUSEMOTION:
        if event.buttons[0]:
          x = self.center[0] - event.rel[0]
          y = self.center[1] - event.rel[1]
          self.center = (x, y)

    return not stop