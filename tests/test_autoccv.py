"""autoccv tests."""

import autoccv


def test_stub() -> None:
    """Stub test to ensure the test suite runs."""
    print(autoccv.__version__)  # noqa: T201


def test__greet() -> None:
    """Test the greet function."""
    assert autoccv.greet("World") == "Hello, World!"


def test__greet_jim() -> None:
    """Test the greet_jim function."""
    assert autoccv.greet_jim() == "Hello, Jim!"
