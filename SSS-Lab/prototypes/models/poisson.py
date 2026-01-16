#!/usr/bin/env python3
"""
Solucion de ∇²u = f(x,y) en un cuadrado unitario con condiciones de 
Dirilecht

Uso:
    python poisson.py --source point --grid 100 --viz heatmap
    python poisson.py --source dipole --viz all
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import argparse


class PoissonSolver2D:

    
    def __init__(self, n=50):
        self.n = n
        self.h = 1.0 / (n - 1)
        self.h2 = self.h * self.h
        self.x = np.linspace(0, 1, n)
        self.y = np.linspace(0, 1, n)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        
    def get_source_term(self, source_type='point'):
        f = np.zeros((self.n, self.n))
        
        if source_type == 'point':
            center = self.n // 2
            f[center, center] = 100.0 / self.h2
            
        elif source_type == 'dipole':
            i1, j1 = self.n // 3, self.n // 2
            i2, j2 = 2 * self.n // 3, self.n // 2
            f[i1, j1] = 100.0 / self.h2
            f[i2, j2] = -100.0 / self.h2
            
        elif source_type == 'ring':
            for i in range(self.n):
                for j in range(self.n):
                    dx = self.x[j] - 0.5
                    dy = self.y[i] - 0.5
                    r = np.sqrt(dx*dx + dy*dy)
                    if 0.2 < r < 0.25:
                        f[i, j] = 50.0 / self.h2
                        
        return f
    
    def solve(self, f, max_iter=5000, tol=1e-6):
        u = np.zeros((self.n, self.n))
        
        for iteration in range(max_iter):
            u_old = u.copy()
            
            # Interior points (Gauss-Seidel)
            for i in range(1, self.n - 1):
                for j in range(1, self.n - 1):
                    u[i, j] = 0.25 * (
                        u[i+1, j] + u[i-1, j] + 
                        u[i, j+1] + u[i, j-1] - 
                        self.h2 * f[i, j]
                    )
            
            # Boundary conditions (Dirichlet: u = 0)
            u[0, :] = 0
            u[-1, :] = 0
            u[:, 0] = 0
            u[:, -1] = 0
            
            # Verificar convergencia 
            if np.max(np.abs(u - u_old)) < tol:
                print(f"{iteration + 1} iterations")
                break
                
        return u
    
    def gradiente(self, u):
        grad_y, grad_x = np.gradient(u, self.h)
        return grad_x, grad_y


class Visualizer
    @staticmethod
    def plot_heatmap(X, Y, u, title="Solution Heatmap"):
        """Plot heatmap of solution"""
        plt.figure(figsize=(8, 7))
        plt.contourf(X, Y, u, levels=50, cmap='RdBu_r')
        plt.colorbar(label='u(x,y)')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(title)
        plt.axis('equal')
        plt.tight_layout()
        
    @staticmethod
    def plot_contours(X, Y, u, title="Equipotential Lines"):
        plt.figure(figsize=(8, 7))
        levels = np.linspace(u.min(), u.max(), 20)
        CS = plt.contour(X, Y, u, levels=levels, cmap='viridis', linewidths=1.5)
        plt.clabel(CS, inline=True, fontsize=8)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(title)
        plt.axis('equal')
        plt.tight_layout()
        
    @staticmethod
    def plot_gradient(X, Y, grad_x, grad_y, title="Gradient Field"):
        plt.figure(figsize=(8, 7))
        
        # Subsample for cleaner visualization
        skip = max(1, len(X) // 20)
        plt.quiver(X[::skip, ::skip], Y[::skip, ::skip], 
                   -grad_x[::skip, ::skip], -grad_y[::skip, ::skip],
                   alpha=0.7, color='blue')
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(title)
        plt.axis('equal')
        plt.tight_layout()
        
    @staticmethod
    def plot_3d_surface(X, Y, u, title="3D Solution Surface"):
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        surf = ax.plot_surface(X, Y, u, cmap='coolwarm', 
                               linewidth=0, antialiased=True)
        
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('u(x,y)')
        ax.set_title(title)
        fig.colorbar(surf, shrink=0.5)
        plt.tight_layout()


def main():
    parser = argparse.ArgumentParser(description='2D Poisson Equation Solver')
    parser.add_argument('--source', type=str, default='point',
                        choices=['point', 'dipole', 'ring'],
                        help='Source configuration')
    parser.add_argument('--grid', type=int, default=50,
                        help='Grid size (n x n)')
    parser.add_argument('--viz', type=str, default='all',
                        choices=['heatmap', 'contour', 'gradient', '3d', 'all'],
                        help='Visualization type')
    parser.add_argument('--save', action='store_true',
                        help='Save plots to files')
    
    args = parser.parse_args()
    
    print(f"Solving 2D Poisson equation...")
    print(f"Grid size: {args.grid}x{args.grid}")
    print(f"Source type: {args.source}")
    
    # Initialize and solve
    solver = PoissonSolver2D(n=args.grid)
    f = solver.get_source_term(args.source)
    u = solver.solve(f)
    grad_x, grad_y = solver.compute_gradient(u)
    
    print(f"Range: [{u.min():.4f}, {u.max():.4f}]")
    
    # Visualize
    viz = Visualizer()
    
    if args.viz == 'all' or args.viz == 'heatmap':
        viz.plot_heatmap(solver.X, solver.Y, u)
        if args.save:
            plt.savefig('heatmap.png', dpi=150)
            
    if args.viz == 'all' or args.viz == 'contour':
        viz.plot_contours(solver.X, solver.Y, u)
        if args.save:
            plt.savefig('contours.png', dpi=150)
            
    if args.viz == 'all' or args.viz == 'gradient':
        viz.plot_gradient(solver.X, solver.Y, grad_x, grad_y)
        if args.save:
            plt.savefig('gradient.png', dpi=150)
            
    if args.viz == 'all' or args.viz == '3d':
        viz.plot_3d_surface(solver.X, solver.Y, u)
        if args.save:
            plt.savefig('surface_3d.png', dpi=150)
    
    plt.show()


if __name__ == '__main__':
    main()
