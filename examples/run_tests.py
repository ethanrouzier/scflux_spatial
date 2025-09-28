#!/usr/bin/env python3
"""
Main test runner for scflux_spatial.

This script provides a convenient way to run all tests with different configurations.
"""

import sys
import subprocess
import argparse
from pathlib import Path


def run_tests(test_type="all", verbose=False, coverage=False):
    """Run tests with specified configuration."""
    
    # Base pytest command
    cmd = ["python", "-m", "pytest"]
    
    # Add verbosity
    if verbose:
        cmd.append("-v")
    else:
        cmd.append("-q")
    
    # Add coverage if requested
    if coverage:
        cmd.extend(["--cov=scflux_spatial", "--cov-report=html", "--cov-report=term"])
    
    # Add test markers based on type
    if test_type == "unit":
        cmd.extend(["-m", "unit"])
    elif test_type == "integration":
        cmd.extend(["-m", "integration"])
    elif test_type == "fast":
        cmd.extend(["-m", "not slow"])
    elif test_type == "slow":
        cmd.extend(["-m", "slow"])
    
    # Add tests directory
    cmd.append("tests/")
    
    print(f"Running command: {' '.join(cmd)}")
    print("-" * 60)
    
    # Run the tests
    try:
        result = subprocess.run(cmd, check=True)
        print("\n✅ All tests passed!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Tests failed with exit code {e.returncode}")
        return False


def run_specific_test(test_file, verbose=False):
    """Run a specific test file."""
    cmd = ["python", "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    else:
        cmd.append("-q")
    
    cmd.append(f"tests/{test_file}")
    
    print(f"Running command: {' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(cmd, check=True)
        print(f"\n✅ Test {test_file} passed!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Test {test_file} failed with exit code {e.returncode}")
        return False


def run_unittest(test_module):
    """Run tests using unittest instead of pytest."""
    cmd = ["python", "-m", "unittest", f"tests.{test_module}", "-v"]
    
    print(f"Running command: {' '.join(cmd)}")
    print("-" * 60)
    
    try:
        result = subprocess.run(cmd, check=True)
        print(f"\n✅ Unittest {test_module} passed!")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n❌ Unittest {test_module} failed with exit code {e.returncode}")
        return False


def main():
    """Main function for test runner."""
    parser = argparse.ArgumentParser(
        description="Test runner for scflux_spatial",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all tests
  python run_tests.py

  # Run only unit tests
  python run_tests.py --type unit

  # Run fast tests only (exclude slow tests)
  python run_tests.py --type fast

  # Run with coverage report
  python run_tests.py --coverage

  # Run specific test file
  python run_tests.py --file test_gpr.py

  # Run using unittest
  python run_tests.py --unittest test_gpr
        """
    )
    
    parser.add_argument(
        "--type", "-t",
        choices=["all", "unit", "integration", "fast", "slow"],
        default="all",
        help="Type of tests to run (default: all)"
    )
    
    parser.add_argument(
        "--file", "-f",
        type=str,
        help="Run specific test file (e.g., test_gpr.py)"
    )
    
    parser.add_argument(
        "--unittest", "-u",
        type=str,
        help="Run specific test module using unittest (e.g., test_gpr)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Verbose output"
    )
    
    parser.add_argument(
        "--coverage", "-c",
        action="store_true",
        help="Generate coverage report"
    )
    
    args = parser.parse_args()
    
    # Check if we're in the right directory
    if not Path("tests").exists():
        print("❌ Error: tests directory not found. Please run from project root.")
        sys.exit(1)
    
    # Run tests based on arguments
    success = False
    
    if args.file:
        success = run_specific_test(args.file, args.verbose)
    elif args.unittest:
        success = run_unittest(args.unittest)
    else:
        success = run_tests(args.type, args.verbose, args.coverage)
    
    if success:
        print("\n🎉 Test run completed successfully!")
        sys.exit(0)
    else:
        print("\n💥 Test run failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
