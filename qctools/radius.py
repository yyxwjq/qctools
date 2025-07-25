# -*- coding: utf-8 -*-
"""
Author: Jiaqi Wang
Date: 2025-07-09 23:45
"""

from typing import List, Tuple, Optional
from collections import deque
import numpy as np
import warnings
import random

class Sphere:
    """Represents a sphere with a center and radius."""
    def __init__(self, center: np.ndarray, radius: float):
        self.center = center
        self.radius = radius
    
    def __str__(self):
        return f"Sphere(center={self.center}, radius={self.radius:.6f})"

class MinimumEnclosingBall:
    """Implements algorithms to find the minimum enclosing ball of a set of points."""
    
    def __init__(self, points: List[np.ndarray], random_seed: Optional[int] = None):
        self.points = np.array(points)
        self.dim = self.points.shape[1] if len(self.points) > 0 else 3
        if random_seed is not None:
            random.seed(random_seed)
            np.random.seed(random_seed)
    
    # ----------------- Helper Methods -----------------
    def _is_inside(self, point: np.ndarray, sphere: Sphere, tolerance: float = 1e-8) -> bool:
        """Checks if a point is inside the sphere with a small tolerance."""
        return np.linalg.norm(point - sphere.center) <= sphere.radius + tolerance
    
    def _validate_sphere(self, sphere: Sphere, tolerance: float = 1e-6) -> bool:
        """Validates if the sphere encloses all points within a tolerance."""
        return np.max(np.linalg.norm(self.points - sphere.center, axis=1)) <= sphere.radius * (1 + tolerance)
    
    def _trivial_sphere(self, support: List[np.ndarray], exact: bool = False) -> Sphere:
        """Computes a sphere for a small number of points (0 to dim+1 points).
        
        Args:
            support: List of points to enclose.
            exact: If True and len(support) == dim+1, compute exact circumsphere using linear system.
        
        Returns:
            Sphere: A sphere enclosing all points in support.
        """
        if not support:
            return Sphere(np.zeros(self.dim), 0)
        if len(support) == 1:
            return Sphere(support[0], 0)
        
        # Use exact method for dim+1 points if requested
        if exact and len(support) == self.dim + 1:
            try:
                # Construct linear system for exact circumsphere
                p1 = support[0]
                A = np.array([p - p1 for p in support[1:]])
                b = 0.5 * np.array([np.sum(p**2) - np.sum(p1**2) for p in support[1:]])
                center = np.linalg.solve(A, b)
                radius = np.linalg.norm(p1 - center)
                return Sphere(center, radius)
            except np.linalg.LinAlgError:
                warnings.warn("Degenerate points detected, falling back to mean method")
        
        # Default: use mean method for 2 to dim+1 points
        center = np.mean(support, axis=0)
        radius = max(np.linalg.norm(p - center) for p in support) if len(support) > 1 else 0
        return Sphere(center, radius)
    
    # ----------------- Welzl Algorithm (Exact) -----------------
    def welzl(self, points: np.ndarray = None, support: List[np.ndarray] = None, exact: bool = False) -> Sphere:
        """Welzl's exact minimum enclosing ball algorithm (recursive).
        
        Args:
            points: The set of points to enclose (default: all points).
            support: Points defining the current sphere's boundary.
            exact: If True, use exact circumsphere computation for dim+1 points.
        
        Returns:
            Sphere: The smallest sphere enclosing all points.
        """
        if points is None:
            points = self.points.copy()
        if support is None:
            support = []
            
        if len(points) == 0 or len(support) == self.dim + 1:
            return self._trivial_sphere(support, exact=exact)
        
        idx = random.randint(0, len(points) - 1)
        p = points[idx]
        points_without_p = np.delete(points, idx, axis=0)
        
        sphere = self.welzl(points_without_p, support, exact=exact)
        return sphere if self._is_inside(p, sphere) else self.welzl(points_without_p, support + [p], exact=exact)
    
    # ----------------- Ritter Algorithm (Approximate) -----------------
    def ritter(self) -> Sphere:
        """Ritter's approximate minimum enclosing ball algorithm."""
        if len(self.points) == 0:
            return Sphere(np.zeros(self.dim), 0)
        
        min_p, max_p = np.argmin(self.points[:, 0]), np.argmax(self.points[:, 0])
        center = (self.points[min_p] + self.points[max_p]) / 2
        radius = np.linalg.norm(self.points[max_p] - center)
        
        for p in self.points:
            dist = np.linalg.norm(p - center)
            if dist > radius:
                radius = (radius + dist) / 2
                center += (p - center) * ((dist - radius) / dist)
        
        return Sphere(center, radius)
    
    # ----------------- Bouncing Bubble Algorithm (Approximate) -----------------
    def bouncing_bubble(self, max_iterations: int = 50) -> Sphere:
        """Bouncing Bubble approximate minimum enclosing ball algorithm."""
        if len(self.points) == 0:
            return Sphere(np.zeros(self.dim), 0)
        
        center, radius = self.points[0].copy(), 0.0
        for _ in range(max_iterations):
            moved = False
            np.random.shuffle(self.points)
            for p in self.points:
                dist = np.linalg.norm(p - center)
                if dist > radius:
                    center += (p - center) * ((dist - radius) / (2 * dist))
                    radius = (radius + dist) / 2
                    moved = True
            if not moved:
                break
        return Sphere(center, radius)
    
    # ----------------- Stable Versions -----------------
    def _run_stable(self, algorithm, max_trials: int = 100, convergence_threshold: float = 0.0005, history_maxlen: int = 20, min_history: int = 10, **kwargs) -> Tuple[Sphere, List[float]]:
        """Runs an algorithm multiple times and selects the best result.
        
        Args:
            algorithm: The algorithm to run (e.g., welzl, bouncing_bubble).
            max_trials: Maximum number of trials to run.
            convergence_threshold: Threshold for convergence check.
            history_maxlen: Maximum length of history deque.
            min_history: Minimum history length for convergence check.
            **kwargs: Additional arguments for the algorithm.
        
        Returns:
            Tuple[Sphere, List[float]]: Best sphere and radius history.
        """
        history, best_sphere, radius_history = deque(maxlen=history_maxlen), None, []
        
        for trial in range(max_trials):
            sphere = algorithm(**kwargs)
            history.append((sphere.radius, sphere.center))
            radius_history.append(sphere.radius)
            
            if best_sphere is None or (sphere.radius < best_sphere.radius and self._validate_sphere(sphere)):
                best_sphere = sphere
            
            if len(history) >= min_history and np.mean(np.abs(np.diff([r for r, _ in list(history)[-min_history:]]))) < np.mean([r for r, _ in list(history)[-min_history:]]) * convergence_threshold:
                print(f"{algorithm.__name__} converged after {trial+1} trials")
                break
        
        # Ensure best_sphere is valid
        if best_sphere is None or not self._validate_sphere(best_sphere):
            warnings.warn("No valid sphere found, returning trivial sphere")
            return self._trivial_sphere(self.points.tolist()), radius_history
        
        return best_sphere, radius_history
    
    def stable_welzl(self, max_trials: int = 100, exact: bool = False) -> Tuple[Sphere, List[float]]:
        """Stable version of Welzl's algorithm."""
        return self._run_stable(self.welzl, max_trials, exact=exact)
    
    def stable_bouncing_bubble(self, max_trials: int = 100, max_iterations: int = 50) -> Tuple[Sphere, List[float]]:
        """Stable version of Bouncing Bubble algorithm."""
        return self._run_stable(self.bouncing_bubble, max_trials, max_iterations=max_iterations)