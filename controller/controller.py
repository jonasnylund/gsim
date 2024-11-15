from collections.abc import Sequence
import time

import numpy as np
import pygame

from gsim.model import py_model  # type: ignore
from gsim.viewer import viewer as viewer_lib


class Controller:
  """Controls the progress of the simulation."""

  def __init__(self, model: py_model.Model, viewer: viewer_lib.Viewer):
    self.model = model
    self.viewer = viewer

    self.running = True

  def update(self) -> bool:
    """Updates the display and polls the events."""
    for event in pygame.event.get():
      if event.type == pygame.QUIT:
        self.running = False
      elif event.type == pygame.MOUSEWHEEL:
        self.viewer.zoom(1 + np.clip(event.y, -2, 2) / 10)
      elif event.type == pygame.MOUSEMOTION:
        if event.buttons[0]:
          self.viewer.move(event.rel)

    return self.running
