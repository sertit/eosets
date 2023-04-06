# """ Testing pairs """
# from sertit import ci
#
# from CI.scripts_utils import data_path
# from eosets.multi_pairs import MultiPairs
#
# ci.reduce_verbosity()
#
#
# def test_pld_pairs():
#     pld_pan = data_path() / "IMG_PHR1A_P_001"
#     pld_ms = data_path() / "IMG_PHR1A_MS_004"
#
#     pairs = MultiPairs(pld_pan, pld_ms)
#
#     ci.assert_val(pairs.has_child, True, "has_child")
#     ci.assert_val(pairs.has_unique_child, True, "has_unique_child")
#     ci.assert_val(pairs.has_children, False, "has_children")
#     ci.assert_val(pairs.same_constellation, True, "same_constellation")
#     ci.assert_val(pairs.same_crs, True, "same_crs")
#     ci.assert_val(pairs.same_sensor_type, True, "same_sensor_type")
#
#
# def test_maxar_pairs():
#     maxar_pan = data_path() / "055670633040_01_P001_PAN"
#     maxar_ms = data_path() / "055670633040_01_P001_MUL"
#
#     pairs = MultiPairs(maxar_pan, maxar_ms)
#
#     ci.assert_val(pairs.has_child, True, "has_child")
#     ci.assert_val(pairs.has_unique_child, True, "has_unique_child")
#     ci.assert_val(pairs.has_children, False, "has_children")
#     ci.assert_val(pairs.same_constellation, True, "same_constellation")
#     ci.assert_val(pairs.same_crs, True, "same_crs")
#     ci.assert_val(pairs.same_sensor_type, True, "same_sensor_type")
#
#
# def test_sv1_pairs():
#     sv1_pan_ms = data_path() / "0032100150001_01"
#
#     pairs = MultiPairs(sv1_pan_ms)
#
#     ci.assert_val(pairs.has_child, False, "has_child")
#     ci.assert_val(pairs.has_unique_child, False, "has_unique_child")
#     ci.assert_val(pairs.has_children, False, "has_children")
#     ci.assert_val(pairs.same_constellation, True, "same_constellation")
#     ci.assert_val(pairs.same_crs, True, "same_crs")
#     ci.assert_val(pairs.same_sensor_type, True, "same_sensor_type")
#
#
# def test_l9_pairs():
#     l9_pan_ms = data_path() / "LC09_L1TP_200030_20220201_20220201_02_T1.tar"
#
#     pairs = MultiPairs(l9_pan_ms)
#
#     ci.assert_val(pairs.has_child, False, "has_child")
#     ci.assert_val(pairs.has_unique_child, False, "has_unique_child")
#     ci.assert_val(pairs.has_children, False, "has_children")
#     ci.assert_val(pairs.same_constellation, True, "same_constellation")
#     ci.assert_val(pairs.same_crs, True, "same_crs")
#     ci.assert_val(pairs.same_sensor_type, True, "same_sensor_type")
