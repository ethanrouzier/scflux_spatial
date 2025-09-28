#!/usr/bin/env python3
"""
Performance tests for scflux_spatial.

This module contains performance benchmarks and scalability tests.
"""

import unittest
import time
import numpy as np
import pytest
from scflux_spatial.gem.gpr import GPRParser


class TestPerformance(unittest.TestCase):
    """Performance tests for core functionality."""
    
    def setUp(self):
        """Set up performance test fixtures."""
        self.gpr_parser = GPRParser()
        
        # Create large datasets for performance testing
        self.large_expression = {}
        for i in range(1000):
            self.large_expression[f'GENE_{i}'] = np.random.uniform(0, 10)
        
        # Complex GPR rules for testing
        self.complex_rules = [
            'GENE_0 and GENE_1 and GENE_2 and GENE_3 and GENE_4',
            '(GENE_0 or GENE_1) and (GENE_2 or GENE_3) and GENE_4',
            '((GENE_0 and GENE_1) or (GENE_2 and GENE_3)) and (GENE_4 or GENE_5)',
            'GENE_0 and (GENE_1 or (GENE_2 and GENE_3)) and (GENE_4 or GENE_5)'
        ]
    
    @pytest.mark.slow
    def test_gpr_evaluation_performance(self):
        """Test GPR evaluation performance with large datasets."""
        # Test with large gene expression dataset
        start_time = time.time()
        
        for rule in self.complex_rules:
            result = self.gpr_parser._evaluate_gpr_rule_with_operators(
                rule, self.large_expression
            )
            self.assertIsInstance(result, (int, float))
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        # Should complete within reasonable time
        self.assertLess(execution_time, 5.0)  # 5 seconds max
        print(f"GPR evaluation time: {execution_time:.3f} seconds")
    
    @pytest.mark.slow
    def test_batch_gpr_evaluation(self):
        """Test batch GPR evaluation performance."""
        # Create multiple rules
        rules = [f'GENE_{i} and GENE_{i+1}' for i in range(0, 100, 2)]
        
        start_time = time.time()
        
        results = {}
        for rule in rules:
            results[rule] = self.gpr_parser._evaluate_gpr_rule_with_operators(
                rule, self.large_expression
            )
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        # Should complete within reasonable time
        self.assertLess(execution_time, 10.0)  # 10 seconds max
        self.assertEqual(len(results), len(rules))
        print(f"Batch GPR evaluation time: {execution_time:.3f} seconds")
    
    def test_memory_usage(self):
        """Test memory usage with large datasets."""
        import psutil
        import os
        
        process = psutil.Process(os.getpid())
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Create large dataset
        large_dataset = {}
        for i in range(10000):
            large_dataset[f'GENE_{i}'] = np.random.uniform(0, 10)
        
        peak_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_increase = peak_memory - initial_memory
        
        # Memory increase should be reasonable
        self.assertLess(memory_increase, 500)  # Less than 500 MB increase
        print(f"Memory increase: {memory_increase:.1f} MB")
    
    @pytest.mark.slow
    def test_scalability_with_dataset_size(self):
        """Test scalability with increasing dataset sizes."""
        sizes = [100, 500, 1000, 2000]
        execution_times = []
        
        for size in sizes:
            # Create dataset of given size
            expression = {}
            for i in range(size):
                expression[f'GENE_{i}'] = np.random.uniform(0, 10)
            
            # Test rule
            rule = ' and '.join([f'GENE_{i}' for i in range(min(10, size))])
            
            start_time = time.time()
            result = self.gpr_parser._evaluate_gpr_rule_with_operators(rule, expression)
            end_time = time.time()
            
            execution_time = end_time - start_time
            execution_times.append(execution_time)
            
            self.assertIsInstance(result, (int, float))
        
        # Execution time should scale reasonably with dataset size
        # (not exponentially)
        max_time = max(execution_times)
        self.assertLess(max_time, 2.0)  # Max 2 seconds
        
        print(f"Execution times by dataset size: {dict(zip(sizes, execution_times))}")
    
    def test_concurrent_evaluation(self):
        """Test concurrent GPR evaluation performance."""
        import threading
        import queue
        
        def evaluate_rule(rule, expression, results_queue):
            """Worker function for concurrent evaluation."""
            start_time = time.time()
            result = self.gpr_parser._evaluate_gpr_rule_with_operators(rule, expression)
            end_time = time.time()
            results_queue.put((rule, result, end_time - start_time))
        
        # Create multiple threads
        threads = []
        results_queue = queue.Queue()
        
        start_time = time.time()
        
        for i, rule in enumerate(self.complex_rules):
            thread = threading.Thread(
                target=evaluate_rule,
                args=(rule, self.large_expression, results_queue)
            )
            threads.append(thread)
            thread.start()
        
        # Wait for all threads to complete
        for thread in threads:
            thread.join()
        
        end_time = time.time()
        total_time = end_time - start_time
        
        # Collect results
        results = []
        while not results_queue.empty():
            results.append(results_queue.get())
        
        # Should complete within reasonable time
        self.assertLess(total_time, 3.0)  # 3 seconds max
        self.assertEqual(len(results), len(self.complex_rules))
        
        print(f"Concurrent evaluation time: {total_time:.3f} seconds")
        print(f"Number of threads: {len(threads)}")


class TestScalability(unittest.TestCase):
    """Scalability tests for different components."""
    
    @pytest.mark.slow
    def test_large_expression_matrix(self):
        """Test performance with large expression matrices."""
        # Simulate large expression matrix (1000 genes x 10000 cells)
        n_genes = 1000
        n_cells = 10000
        
        # Create expression matrix
        start_time = time.time()
        expression_matrix = np.random.poisson(5, (n_cells, n_genes))
        end_time = time.time()
        
        creation_time = end_time - start_time
        self.assertLess(creation_time, 5.0)  # 5 seconds max
        
        # Test operations on large matrix
        start_time = time.time()
        
        # Calculate mean expression per gene
        mean_expression = np.mean(expression_matrix, axis=0)
        
        # Calculate variance per gene
        var_expression = np.var(expression_matrix, axis=0)
        
        # Find highly variable genes
        high_var_genes = np.where(var_expression > np.percentile(var_expression, 90))[0]
        
        end_time = time.time()
        operation_time = end_time - start_time
        
        self.assertLess(operation_time, 2.0)  # 2 seconds max
        self.assertEqual(len(mean_expression), n_genes)
        self.assertEqual(len(high_var_genes), n_genes // 10)  # Top 10%
        
        print(f"Matrix creation time: {creation_time:.3f} seconds")
        print(f"Matrix operations time: {operation_time:.3f} seconds")
    
    def test_memory_efficient_operations(self):
        """Test memory-efficient operations on large datasets."""
        # Create large dataset
        large_data = np.random.uniform(0, 10, (1000, 1000))
        
        # Test memory-efficient operations
        start_time = time.time()
        
        # Use memory-efficient operations
        result = np.sum(large_data, axis=1)  # Sum along rows
        result = np.mean(large_data, axis=0)  # Mean along columns
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        self.assertLess(execution_time, 1.0)  # 1 second max
        self.assertEqual(len(result), 1000)
        
        print(f"Memory-efficient operations time: {execution_time:.3f} seconds")


if __name__ == '__main__':
    unittest.main(verbosity=2)
