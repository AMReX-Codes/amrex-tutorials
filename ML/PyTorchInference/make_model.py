#!/usr/bin/env python3
"""Write model.pt, the TorchScript model used by main.cpp.

The model is one dense layer, y = tanh(W x + b), with fixed weights so that
main.cpp can check the result. A trained network is saved the same way.
"""

import torch


class Model(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.linear = torch.nn.Linear(2, 2)
        with torch.no_grad():
            # Must match the weights in main.cpp
            self.linear.weight.copy_(torch.tensor([[1.0, -0.5], [0.25, 2.0]]))
            self.linear.bias.copy_(torch.tensor([0.5, -0.25]))

    def forward(self, x):
        # x has shape (number of cells, 2)
        return torch.tanh(self.linear(x))


if __name__ == "__main__":
    # Record the operations for an example batch; any batch size works later.
    model = torch.jit.trace(Model().eval(), torch.zeros(8, 2))
    model.save("model.pt")
    print("Wrote model.pt")
