import json
import tempfile
import unittest
from pathlib import Path

import surface_stress


class SurfaceStressTests(unittest.TestCase):
    def test_pressure_conversion_and_sample_folder(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            info = root / "sample.info"
            info.write_text(json.dumps({
                "geometry": "film",
                "case_name": "sample",
                "simulation_template": {"surface_measurement": {"enabled": True}},
            }), encoding="utf-8")
            stress = root / "stress.sample.surface.dat"
            stress.write_text(
                "# time_fs temp_K pe_kcal_per_mol pxx_atm pyy_atm pzz_atm lx_A ly_A lz_A\n"
                "5000 300 -10 0 0 2 20 20 100\n"
                "10000 300 -10 0 0 2 20 20 100\n",
                encoding="utf-8",
            )
            result = surface_stress.analyze(info, stress, root / "sample", 0.005)
            self.assertAlmostEqual(result["mean_sample_mN_per_m"], 1.01325)
            self.assertEqual(result["samples"], 2)
            self.assertTrue((root / "sample" / "surface_series.sample.dat").is_file())
            self.assertTrue((root / "sample" / "surface_blocks.sample.dat").is_file())
            self.assertTrue((root / "sample" / "surface_summary.sample.json").is_file())

    def test_non_film_info_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            info = Path(temporary) / "bad.info"
            info.write_text('{"geometry":"bulk","case_name":"bad"}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "film"):
                surface_stress.load_info(info)


if __name__ == "__main__":
    unittest.main()
