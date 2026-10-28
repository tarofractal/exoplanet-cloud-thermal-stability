"""Analytic and synthetic checks; no publication inputs are used."""
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / 'scripts' / f'{name}.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


qabs = load_script('qabs_loggrid')
thermal = load_script('thermal_decoupling')


class NumericalChecks(unittest.TestCase):
    def synthetic_spectrum(self, directory, n, k, radius):
        # Deliberately artificial constant-index material for analytic checks.
        np.savetxt(Path(directory) / 'Synthetic.dat', [[1.0, n, k], [2.0, n, k]])
        return qabs.compute_qabs_on_log_grid(
            'Synthetic', directory, [radius], R=128,
            save_dir=directory, make_plot=False,
        )[0]

    def test_small_absorbing_sphere_matches_rayleigh_limit(self):
        radius = 0.001
        with tempfile.TemporaryDirectory() as directory:
            spectrum = self.synthetic_spectrum(directory, 1.5, 0.1, radius)
            wavelength = spectrum['Wavelength (μm)'].to_numpy()
            m = complex(1.5, 0.1)
            expected = 4 * (2 * np.pi * radius / wavelength) * ((m*m-1)/(m*m+2)).imag
            np.testing.assert_allclose(spectrum['Qabs_0.001μm'], expected, rtol=1e-3, atol=1e-12)

    def test_nonabsorbing_sphere_has_zero_absorption(self):
        with tempfile.TemporaryDirectory() as directory:
            spectrum = self.synthetic_spectrum(directory, 1.5, 0.0, 0.1)
            np.testing.assert_allclose(spectrum['Qabs_0.1μm'], 0.0, atol=1e-10)

    def test_constant_efficiency_has_constant_planck_mean(self):
        wavelength = np.geomspace(0.1, 10000, 2000)
        for temperature in (1000.0, 6000.0):
            mean = qabs.integrate_qabs(wavelength, np.full_like(wavelength, 0.37), temperature)
            self.assertAlmostEqual(mean['planck_weighted'], 0.37, places=12)

    def test_qabs_output_can_be_read_by_decoupling_code(self):
        with tempfile.TemporaryDirectory() as directory:
            spectrum = self.synthetic_spectrum(directory, 1.5, 0.1, 0.1)
            wavelength_m, efficiencies = thermal.load_qabs_table_exact(
                Path(directory) / 'Synthetic_logR128_Qabs.csv'
            )
            np.testing.assert_allclose(wavelength_m, spectrum['Wavelength (μm)'] * 1e-6)
            np.testing.assert_allclose(efficiencies[0.1], spectrum['Qabs_0.1μm'])

    def test_grey_particle_has_planetary_equilibrium_temperature(self):
        wavelength_m = np.geomspace(1e-9, 1e-1, 2500)
        efficiencies = np.ones_like(wavelength_m)
        stellar_temperature, equilibrium_temperature = 6000.0, 1000.0
        heating = 4 * (equilibrium_temperature/stellar_temperature)**4 * trapezoid(
            thermal.planck_lambda(wavelength_m, stellar_temperature), wavelength_m
        )
        emission = thermal.EmissionInterpolator(wavelength_m, efficiencies, thermal.EMISSION_T_GRID)
        particle_temperature = thermal.solve_tp(heating, equilibrium_temperature, 0, emission, collision=False)
        self.assertAlmostEqual(particle_temperature, equilibrium_temperature, delta=0.2)

    def test_collisions_drive_particle_temperature_toward_gas(self):
        wavelength_m = np.geomspace(1e-9, 1e-1, 2500)
        emission = thermal.EmissionInterpolator(wavelength_m, np.ones_like(wavelength_m), thermal.EMISSION_T_GRID)
        heating = float(emission(1000.0))
        gas_temperature = 1500.0
        temperatures = [thermal.solve_tp(heating, gas_temperature, p, emission) for p in (0, 1e-6, 100)]
        self.assertTrue(all(1000-0.01 <= t <= gas_temperature for t in temperatures))
        self.assertTrue(temperatures[0] < temperatures[1] < temperatures[2])
        self.assertLess(gas_temperature-temperatures[-1], 1.0)

    def test_missing_input_is_a_failed_command(self):
        with tempfile.TemporaryDirectory() as directory:
            result = subprocess.run([
                sys.executable, str(ROOT / 'scripts/qabs_loggrid.py'),
                '--materials', 'Synthetic', '--data_dir', directory,
                '--save_dir', directory, '--no_plot',
            ], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('No such file', result.stdout)
            self.assertFalse(list(Path(directory).glob('*Qabs.csv')))


if __name__ == '__main__':
    unittest.main()
