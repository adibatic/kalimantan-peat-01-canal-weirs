"""Closed-form relations: each one checked against its own definition."""
import numpy as np
import pytest

from canalweir.hydraulics import (G, WEIR_COEFF, critical_depth, froude, hydraulic_radius,
                                  manning_unit_discharge, normal_depth, weir_head,
                                  weir_unit_discharge)


def test_hydraulic_radius_of_a_rectangle():
    # 3 m wide, 1 m deep: area 3 m^2, wetted perimeter 5 m
    assert hydraulic_radius(1.0, 3.0) == pytest.approx(0.6, rel=1e-15)


def test_wide_canal_limit():
    assert hydraulic_radius(0.7, 1e9) == pytest.approx(0.7, rel=1e-8)
    assert hydraulic_radius(0.7) == 0.7


@pytest.mark.parametrize("width", [None, 3.0])
def test_normal_depth_inverts_manning(width):
    q, slope, n = 0.2, 2e-4, 0.03
    h = normal_depth(q, slope, n, width)
    assert manning_unit_discharge(h, slope, n, width) == pytest.approx(q, rel=1e-12)


def test_critical_depth_has_froude_one():
    q = 0.35
    assert froude(critical_depth(q), q) == pytest.approx(1.0, rel=1e-13)


def test_weir_coefficient_from_critical_flow_on_the_crest():
    # critical depth on the crest is 2/3 of the head: q = sqrt(G (2H/3)^3)
    head = 0.27
    q = np.sqrt(G * (2.0 * head / 3.0) ** 3)
    assert weir_unit_discharge(10.0 + head, 10.0) == pytest.approx(q, rel=1e-14)
    assert WEIR_COEFF == pytest.approx(np.sqrt(8.0 * G / 27.0), rel=1e-15)


def test_weir_head_inverts_weir_law():
    q = 0.123
    assert weir_unit_discharge(weir_head(q), 0.0) == pytest.approx(q, rel=1e-13)


def test_no_flow_below_the_crest():
    assert weir_unit_discharge(0.99, 1.0) == 0.0
    assert weir_unit_discharge(1.0, 1.0) == 0.0