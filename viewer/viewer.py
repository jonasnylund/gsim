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
      flags=pygame.RESIZABLE,
    )
    self.scale: float = 1
    self.center: tuple[float, float] = (500, 250)
    self.rendering_time: float = 0

  def __del__(self):
    pygame.quit()

  def show(self, particles: Sequence[py_model.Particle]) -> None:
    """Renders all the given particles to the display."""
    start_time = time.time()

    self.display.fill((0, 0, 0))
    for particle in particles:
      x, y = particle.position[:2]
      size = np.pow(particle.mass, 1 / 3) * 2 * self.scale
      opacity = min(size, 2) / 2 * 255
      center=(
          x * self.scale + self.center[0],
          y * self.scale + self.center[1],
      )
      draw_circle_alpha(
        self.display,
        color=(200, 200, 255, int(opacity)),
        center=center,
        radius=max(1, size),
      )
    pygame.display.flip()

    self.rendering_time += time.time() - start_time

  def zoom(self, factor) -> None:
    """Changes the scale of the view."""
    self.scale *= factor

  def move(self, delta: tuple[float, float]) -> None:
    """Move the center point in the view."""
    self.center = (
      self.center[0] + delta[0],
      self.center[1] + delta[1],
    )
