"""Temperature-range parsing and midpoint-based green scoring."""

from __future__ import annotations


def parse_temperature_range(temp: str | None) -> tuple[float | None, float | None]:
    """ Parses a temperature range from a string in the format "(lo, hi)". Returns a tuple of (lo, hi) as floats, or (None, None) if the input is invalid. """
    if temp is None or not isinstance(temp, str):
        return None, None
    temp = temp.strip()
    if not temp:
        return None, None
    try:
        lo_str, hi_str = [item.strip() for item in temp[1:-1].split(",")]
        return float(lo_str), float(hi_str)
    except Exception:
        return None, None


def temperature_score_from_range(temp: str | None) -> float | None:
    """ Calculates a temperature score based on a temperature range. Returns the score as a float, or None if the input is invalid. """
    lo, hi = parse_temperature_range(temp)
    if lo is None or hi is None:
        return None

    midpoint = (lo + hi) / 2.0
    if 10 <= midpoint <= 30:
        return 1.00
    if 0 <= midpoint <= 60:
        return 0.75
    if -10 <= midpoint <= 100:
        return 0.50
    if -30 <= midpoint <= 150:
        return 0.25
    return 0.10
