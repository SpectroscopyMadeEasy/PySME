"""No EOS or stellar synthesis; requires CFLD_LIBRARY pointing to isolated build."""
import importlib.util
import os
from pathlib import Path
import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('cfld_api',Path(__file__).resolve().parents[2]/'src/pysme/line_depression.py')
api=importlib.util.module_from_spec(spec);spec.loader.exec_module(api)
@pytest.fixture
def native():
    path=os.environ.get('CFLD_LIBRARY')
    if not path:pytest.skip('Set CFLD_LIBRARY to isolated diagnostic library')
    return api.NativeLineDepression(path)

@pytest.mark.parametrize('tangent',[False,True])
def test_closure(native,tangent):
    m=np.geomspace(.001,20,64);b=1+.5*m;k=1+.2*np.sqrt(m)
    r=native.prescribed_path(m,k,k+.7*np.exp(-m/10),b+.05*m*m,b,b,tangent=tangent,quadrature=32)
    np.testing.assert_allclose(r.total.sum()+r.boundary,r.metadata[0]-r.metadata[1],rtol=1e-6)
    for ch,bd,val in [(2,2,0),(3,3,1)]:
        np.testing.assert_allclose(r.channels[ch].sum()+r.metadata[bd],r.metadata[val],rtol=1e-14)
    np.testing.assert_allclose(r.channels[1],r.channels[4]+r.channels[5],rtol=1e-12,atol=1e-10)
    x=np.log10(m)
    np.testing.assert_allclose(r.density(x)*np.diff(x),r.total)
    np.testing.assert_allclose(r.density(-x)*abs(np.diff(-x)),r.total)


def test_no_line(native):
    m=np.geomspace(.001,20,56);b=1+m;k=np.ones(56)
    r=native.prescribed_path(m,k,k,b,b,b)
    np.testing.assert_array_equal(r.total,0)
    assert r.boundary==0

@pytest.mark.parametrize('bad',[np.nan,np.inf,-1.,0.])
def test_bad_opacity(native,bad):
    m=np.arange(3.);k=np.ones(3);b=np.ones(3);k[1]=bad
    with pytest.raises(ValueError):native.prescribed_path(m,k,b,b,b,b)


def test_unprepared_state(native):
    with pytest.raises(RuntimeError,match='status 10'):
        native.spherical_state([5000.],[1.],56)


def test_global_weights_not_chunk_normalised():
    r=api.LineDepressionResult(np.ones((2,3,6,4)),np.ones((2,3,9)))
    whole=r.disc_flux([.1,.3,.6])
    parts=[api.LineDepressionResult(r.channels[:,i:i+1],r.metadata[:,i:i+1]).disc_flux([w]) for i,w in enumerate([.1,.3,.6])]
    for k in [0,1]:np.testing.assert_allclose(sum(x[k] for x in parts),whole[k])

@pytest.mark.parametrize('tangent', [False, True])
def test_signed_source_and_repeated_calls(native, tangent):
    """A changed intervening call must not leak its trace into the next call."""
    m = np.geomspace(.001, 20, 56)
    k = np.ones_like(m)
    b = 1 + .5*m
    args = (m, k, 2*k, b, b, b)
    first = native.prescribed_path(*args, tangent=tangent)
    altered = native.prescribed_path(m, k, 2*k, .5*b, b, b,
                                    tangent=not tangent)
    assert np.any(altered.channels[4] < 0)
    np.testing.assert_allclose(altered.total.sum() + altered.boundary,
                               altered.metadata[0] - altered.metadata[1], rtol=1e-6)
    repeated = native.prescribed_path(*args, tangent=tangent)
    np.testing.assert_array_equal(first.channels, repeated.channels)
    np.testing.assert_array_equal(first.metadata, repeated.metadata)

@pytest.mark.parametrize('mass', [[0., 1.], [0., 1., 1.], [2., 1., 0.]])
def test_invalid_path_grid(native, mass):
    k = np.ones(len(mass))
    with pytest.raises(ValueError):
        native.prescribed_path(mass, k, k, k, k, k)

@pytest.mark.parametrize('quadrature', [3, 129])
def test_invalid_quadrature(native, quadrature):
    k = np.ones(3)
    with pytest.raises(ValueError):
        native.prescribed_path([0., 1., 2.], k, k, k, k, k,
                               quadrature=quadrature)

@pytest.mark.parametrize('edges', [[0., 1., 1.], [0., 2., 1.], [0., 1.]])
def test_invalid_density_coordinates(edges):
    result = api.LineDepressionResult(np.ones((6, 2)), np.ones(7))
    with pytest.raises(ValueError):
        result.density(edges)

@pytest.mark.parametrize('weights', [[-.1, 1.1], [1.], [np.nan, 1.]])
def test_invalid_disc_weights(weights):
    result = api.LineDepressionResult(np.ones((1, 2, 6, 3)), np.ones((1, 2, 9)))
    with pytest.raises(ValueError):
        result.disc_flux(weights)
