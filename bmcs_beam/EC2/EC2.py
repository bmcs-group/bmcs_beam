# # Import Standard Libraries
# from abc import ABCMeta, abstractmethod
# import logging
# import scipy as np
#
# # Import Local Libraries
# import Util_ACI as ACI
# import Util_EC2 as EC2
#
#
# # ===========================================================================
# #   Reinforced Beam
# # ===========================================================================
#
# class ReinforcedBeam(metaclass=ABCMeta):
#     """ Abstract Reinforced Beam Class
#
#     Static Members:
#         gamma_c = 1.5
#         gamma_S = 1.15
#         niter = 20
#     Instance Members:
#
#     Static Methods:
#         plotArrow(self, ax, x, y, dx, dy, color)
#         plotStress(self, ax, fcd, h_top, h_bot)
#
#     """
#
#     gamma_c = 1.5
#     gamma_S = 1.15
#     niter = 20
#
#     H = d_S = As = fck = fyk = Es = k = 0
#     units = "MPa"
#
#     def __init__(self, name=None):
#         self.ta = self.H / 500
#         self.name = name
#
#     def plotArrow(self, ax, x, y, dx, dy, color):
#         if x == 0:
#             ax.arrow(x, y, dx - self.ta, dy, fc=color, ec=color,
#                      width=10 * self.ta, head_width=10 * self.ta, head_length=self.ta)
#         else:
#             ax.arrow(x, y, dx + self.ta, dy, fc=color, ec=color,
#                      width=self.ta, head_width=10 * self.ta, head_length=self.ta)
#
#     @staticmethod
#     def plotStress(ax, fcd, h_top, h_bot):
#         X0 = [-fcd, 0]
#         Y0 = [h_top, h_top]
#         Y1 = [h_bot, h_bot]
#         ax.fill_between(X0, Y0, Y1, color='b', alpha=0.3)
#
#     @staticmethod
#     def neutralAxis(h, bb, bt):
#         """ Input:	h = trapezoidal height
#                     bb = bottom width
#                     bt = top width
#             Output:	yc = bottom to neutral axis """
#         return h * (2 * bt + bb) / 3 / (bb + bt)
#
#     @staticmethod
#     def momentInertia(h, bb, bt):
#         """ Input:	h = trapezoidal height
#                     bb = bottom width
#                     bt = top width
#             Output:	I = moment of inertia about neutral axis """
#         return h ** 3 * (bt ** 2 + 4 * bt * bb + bb ** 2) / 36 / (bt + bb)
#
#     @staticmethod
#     def beamFactory():
#         pass
#
#
# # ===========================================================================
# #   Rectangular Beam
# # ===========================================================================
#
#
# class RectangularBeam(ReinforcedBeam):
#     # b = 558.8
#     # h = 609.6
#     # d = 546.1
#     # As = 3870
#     # fck = 27.6
#     # fyk = 414
#     # rho = 25
#     # Es = 200000
#     # units = "MPa"
#     b = 22
#     h = 24
#     d = 21.5
#     As = 6
#     fck = 4000
#     fyk = 60000
#     rho = 145
#     Es = 29000000
#     units = "psi"
#
#     # ---------------------------------------------------------------------------
#     #   ACI Equations
#     # ---------------------------------------------------------------------------
#
#     def ACI_cracking_moment(self):
#         logging.debug("Uncracked moment capacity per ACI")
#         Ec = ACI.elastic_modulus(self.fck, self.rho, self.units)
#         fr = ACI.tensile_strength(self.fck, self.units)
#         n = self.Es / Ec
#         logging.debug("    fr = {:6.2f}, Es/Ec = {:3.2f}".format(fr, n))
#
#         Ac = self.b * self.h
#         As = (n - 1) * self.As
#         logging.debug("    b*h = {:6.2f}, (n-1)*As = {:6.2f}".format(Ac, As))
#
#         y_bott = (Ac * self.h / 2 + As * (self.h - self.d)) / (Ac + As)
#         I_uncr = self.b * self.h ** 3 / 12 \
#                  + Ac * (self.h / 2 - y_bott) ** 2 \
#                  + As * (y_bott - (self.h - self.d)) ** 2
#         logging.debug("    Botton to NA: y_bott = {:6.2f}".format(y_bott))
#         logging.debug("    Moment of inertia = {:10.2f} ".format(I_uncr))
#
#         Mcr = fr * I_uncr / y_bott
#         logging.info("    Mcr = {:5.2f}".format(Mcr))
#         return Mcr
#
#     def ACI_elastic_moment(self):
#         Ec = ACI.elastic_modulus(self.fck, self.rho, self.units)
#         fc = 0.5 * self.fck
#         n = self.Es / Ec
#         As = n * self.As
#         fs = self.fyk / n
#         logging.debug("    n = {:3.2f}, n*As = {:6.2f}".format(n, As))
#
#         kd = (-As + np.sqrt(As ** 2 + 2 * self.b * As * self.d)) / (self.b)
#         Icr = self.b * kd ** 3 / 12 + self.b * kd * (kd / 2) ** 2 + As * (self.d - kd) ** 2
#         logging.debug("    Top to NA: kd = {:6.2f}".format(kd))
#         logging.debug("    Moment of inertia = {:10.2f}".format(Icr))
#
#         Mel = min(fc * Icr / kd, fs * Icr / (self.d - kd))
#         logging.info("    Mel = {:5.2f}".format(Mel))
#         return Mel
#
#     def ACI_design_moment(self):
#         beta1 = ACI.beta(self.fck, self.units)
#         logging.debug("    beta1 = {:3.2f}".format(beta1))
#         c = (self.As * self.fyk) / (0.85 * self.fck * self.b) / beta1
#         logging.debug("    c = {:6.2f}".format(c))
#         phi = ACI.ductility_requirement(c, self.d, type="beam")
#         logging.debug("    phi = {:3.2f}".format(phi))
#         ACI.steel_ratio(self.As, self.fck, self.fyk,
#                         self.b, self.d, self.units)
#         MRd = phi * self.As * self.fyk * (self.d - beta1 * c / 2)
#         logging.info("    MRd = {:5.2f}".format(MRd))
#         return MRd
#
#     # ---------------------------------------------------------------------------
#     #   EC2 Equations
#     # ---------------------------------------------------------------------------
#
#     def EC2_cracking_moment(self):
#         logging.debug("Uncracked moment capacity per EC2")
#         Ec = EC2.elastic_modulus(self.fck, self.units)
#         fr = EC2.flex_tensile_strength(self.fck, self.h, self.units)
#         n = self.Es / Ec
#         logging.debug("    fr = {:6.2f}, Es/Ec = {:3.2f}".format(fr, n))
#
#         Ac = self.b * self.h
#         As = (n - 1) * self.As
#         logging.debug("    b*h = {:6.2f}, (n-1)*As = {:6.2f}".format(Ac, As))
#
#         y_bott = (Ac * self.h / 2 + As * (self.h - self.d)) / (Ac + As)
#         I_uncr = self.b * self.h ** 3 / 12 \
#                  + Ac * (self.h / 2 - y_bott) ** 2 \
#                  + As * (y_bott - (self.h - self.d)) ** 2
#         logging.debug("    Botton to NA: y_bott = {:6.2f}".format(y_bott))
#         logging.debug("    Moment of inertia = {:10.2f} ".format(I_uncr))
#
#         Mcr = fr * I_uncr / y_bott
#         logging.info("    Mcr = {:5.2f}".format(Mcr))
#         return Mcr
#
#     def EC2_elastic_moment(self):
#         Ec = EC2.elastic_modulus(self.fck, self.units)
#         fc = 0.5 * self.fck
#         n = self.Es / Ec
#         As = n * self.As
#         fs = self.fyk / n
#         logging.debug("    n = {:3.2f}, n*As = {:6.2f}".format(n, As))
#
#         kd = (-As + np.sqrt(As ** 2 + 2 * self.b * As * self.d)) / (self.b)
#         Icr = self.b * kd ** 3 / 12 + self.b * kd * (kd / 2) ** 2 + As * (self.d - kd) ** 2
#         logging.debug("    Top to NA: kd = {:6.2f}".format(kd))
#         logging.debug("    Moment of inertia = {:10.2f}".format(Icr))
#
#         Mel = min(fc * Icr / kd, fs * Icr / (self.d - kd))
#         logging.info("    Mel = {:5.2f}".format(Mel))
#         return Mel
#
#     def EC2_design_moment(self):
#         [alpha, beta] = EC2.alpha_beta(self.fck, self.units)
#         logging.debug("    a = {:3.2f}, b = {:3.2f}".format(alpha, beta))
#         fcd = self.fck / self.gamma_c
#         fyd = self.fyk / self.gamma_S
#         logging.debug("    fcd = {:6.2f}, fyd = {:6.2f}".format(fcd, fyd))
#         Xu = self.As * fyd / (alpha * self.b * fcd)
#         Xu_max = EC2.ductility_requirement(
#             Xu, self.d, self.fck, fyd, self.units)
#
#         EC2.steel_ratio(self.As, self.fck, self.fyk, self.b,
#                         self.d, self.h, Xu_max, self.units)
#
#         MRd = self.As * fyd * (self.d - beta * Xu)
#         logging.info("    MRd = {:5.2f}".format(MRd))
#         return MRd
#
#
# class DoublyReinforcedBeam():
#     a = 2
#     h = 609.6
#     d_1 = 20
#     As_1 = 10
#     d_2 = 546.1
#     As_2 = 40
#     b = 558.8
#     As = 2.8575 ** 2 * np.pi / 4
#     fck = 27.6
#     fyk = 414
#     rho = 25
#     Es = 200000
#     units = "MPa"
#
#
# if __name__ == "__main__":
#     logging.basicConfig(level=logging.DEBUG, format='%(message)s')
#
#     # print("\n Ec for C30/37: \n")
#     # print(ACI.elastic_modulus(30, rho=25, units="MPa"))
#     # print(EC2.elastic_modulus(30))
#     # print("\n Tensile strength for C30/37: \n")
#     # print(ACI.tensile_strength(30, units="MPa"))
#     # print(EC2.tensile_strength(30))
#     # print("\n Tensile strength for C55/67: \n")
#     # print(ACI.tensile_strength(55, units="MPa"))
#     # print(EC2.tensile_strength(55))
#     # print("\n Ultimate strain for C30/37 \n")
#     # print(ACI.ultimate_strain(30, units="MPa"))
#     # print(EC2.ultimate_strain(30))
#     # print("\n Ultimate strain for C55/67 \n")
#     # print(ACI.ultimate_strain(55, units="MPa"))
#     # print(EC2.ultimate_strain(55))
#     # print("\n alpha & Beta factors: \n")
#     # print(EC2.alpha_beta(10))
#     # print(EC2.alpha_beta(60))
#     # print(EC2.alpha_beta(80))
#     # print("\n Lambda & Eta factors: \n")
#     # print(EC2.lambda_eta(10))
#     # print(EC2.lambda_eta(60))
#     # print(EC2.lambda_eta(80))
#
#     beam = RectangularBeam()
#     # print('\n Cracking Moment: \n')
#     # MRd = beam.ACI.cracking_moment()
#     # print(MRd/12000)
#     # MRd = beam.EC2.cracking_moment()
#     # print(MRd/12000)
#
#     print('\n Elastic Moment: \n')
#     MRd = beam.ACI_elastic_moment()
#     print(MRd / 12000)
#     MRd = beam.EC2_elastic_moment()
#     print(MRd / 12000)
#
#     print('\n Design Moment: \n')
#     MRd = beam.ACI_design_moment()
#     print(MRd / 12000)
#     MRd = beam.EC2_design_moment()
#     print(MRd / 12000)
#
#
#
#     ------------------------------------
#
#     # Import Standard Libraries
#     import logging
#     import scipy as np
#
#     # Import Local Libraries
#     from Utilities import *
#
#
#     # ===========================================================================
#     #   EC2 Equations - Material properties
#     # ===========================================================================
#
#     def elastic_modulus(fck, units="MPa"):
#         """ Input:  fck = char. comp. strength of concrete
#                     units = "MPa" or "psi" (default = "MPa")
#             Output: Ec = mean elastic modulus of concrete """
#         fck = convert_2_MPa(fck, units)
#         fcm = fck + 8
#         Ec = 22000 * (fcm / 10) ** 0.3
#         return Ec if units == "MPa" else convert_2_psi(Ec, "MPa")
#
#
#     def tensile_strength(fck, units="MPa"):
#         """ Input:  fck = char. comp. strength of concrete
#                     units = "MPa" or "psi" (default = "MPa")
#             Output: fctm = mean tensile strength of concrete """
#         fck = convert_2_MPa(fck, units)
#         fcm = fck + 8
#         fctm = 0.3 * fck ** (2 / 3) if fck <= 50 else 2.12 * np.log(1 + fcm / 10)
#         return fctm if units == "MPa" else convert_2_psi(fctm, "MPa")
#
#
#     def flex_tensile_strength(fck, h, units="MPa"):
#         """ Input:  fck = char. comp. strength of concrete
#                     h = height of reinforced concrete beam
#                     units = "MPa" or "psi" (default = "MPa")
#             Output: fctm,fl = mean tensile strength for flexure """
#         fck = convert_2_MPa(fck, units)
#         fctm = tensile_strength(fck)
#         h = convert_2_mm(h, units)
#         fctm = min((1.6 - h / 1000) * fctm, fctm)
#         return fctm if units == "MPa" else convert_2_psi(fctm, "MPa")
#
#
#     def ultimate_strain(fck, units="MPa"):
#         """ Input:  fck = char. comp. strength of concrete
#                     units = "MPa" or "psi" (default = "MPa")
#             Output: ecu3 = ultimate tensile strain """
#         fck = convert_2_MPa(fck, units)
#         ecu3 = 2.6 + 35 * ((90 - fck) / 100) ** 4
#         return min(ecu3, 3.5) / 1000
#
#
#     # ===========================================================================
#     #   EC2 Equations - Parameters
#     # ===========================================================================
#
#     def alpha_beta(fck, units="MPa"):
#         """ Input:  fck = char. comp. strength of concrete
#                     units = "MPa" or "psi" (default = "MPa")
#             Output: alpha = factor for bilinear stress block
#                     beta = (dist. from comp. to Nc)/Xu """
#         fck = convert_2_MPa(fck, units)
#         alpha = np.ceil((9E-05 * fck ** 2 - 0.0177 * fck + 1.4032) * 100) / 100
#         beta = np.ceil((4E-05 * fck ** 2 - 0.0071 * fck + 0.634) * 100) / 100
#         return [min(alpha, 0.75), min(beta, 0.39)]
#
#
#     def lambda_eta(fck, units="MPa"):
#         """ Input:  fck = char. comp. strength of concrete
#                     units = "MPa" or "psi" (default = "MPa")
#             Output: la = (height of compressive zone)/Xu
#                     eta = factor for "Whitney" stress block """
#         fck = convert_2_MPa(fck, units)
#         la = min(0.8 - (fck - 50) / 400, 0.8)
#         eta = min(1 - (fck - 50) / 200, 1.0)
#         return [la, eta]
#
#
#     # ===========================================================================
#     #   EC2 Equations - Maximum reinforcement (Ductility)
#     # ===========================================================================
#
#     def ductility_requirement(Xu, d, fck, fyd, units="MPa"):
#         """ Input:  Xu = dist. from comp. to neutral axis
#                     d = dist. from comp. to reinforcement
#                     fck = char. comp. strength of concrete
#                     fyd = design steel yield stress
#                     units = "MPa" or "psi" (default = "MPa")
#             Output: Xu_max = Max. dist. to neutral axis """
#         [fck, fyd] = convert_2_MPa(np.array([fck, fyd]), units)
#
#         ecu = ultimate_strain(fck)  # units="MPa"
#         Xu_max = min(ecu * 10 ** 6 / (ecu * 10 ** 6 + 7 * fyd), 0.535) * d
#         if Xu < Xu_max:
#             logging.info(
#                 "    Xu = {:6.2f} < Xu_max = {:6.2f}. OK".format(Xu, Xu_max))
#         else:
#             logging.info(
#                 "    Xu = {:6.2f} > Xu_max = {:6.2f}. Not OK".format(Xu, Xu_max))
#         return Xu_max
#
#
#     # ===========================================================================
#     #   EC2 Equations - Minimum reinforcement (Md > Mcr)
#     # ===========================================================================
#
#     def steel_ratio(As, fck, fyk, b, d, h, Xu, units="MPa"):
#         """ Input:  As = area of reinforcement steel
#                     fck = char. comp. strength of concrete
#                     fyk = char. yield stress of reinforcement
#                     b = width of beam portion in compression
#                     d = dist. from comp. to reinforcement
#                     h = height of reinforced concrete beam
#                     Xu = maximum dist. to neutral axis
#                     units = "MPa" or "psi" (default = "MPa")
#             Output: A_min = minimum reinforcement area
#                     A_max = maximum reinforcement area """
#         [fck, fyk] = convert_2_MPa(np.array([fck, fyk]), units)
#         [b, d, h, Xu] = convert_2_mm(np.array([b, d, h, Xu]), units)
#         As = convert_2_mm2(As, units)
#
#         fctm = flex_tensile_strength(fck, h)  # units="MPa"
#         A_min = max((0.26 * fctm / fyk), 0.0013) * (b * d)
#
#         fcd = fck / 1.5
#         fyd = fyk / 1.15
#
#         alpha = alpha_beta(fck)[0]  # units="MPa"
#         A_max = min(alpha * (fcd / fyd) * b * Xu, 0.4 * b * d)
#
#         compare_steel_area(As, A_min, A_max)
#
#         return [A_min, A_max] if units == "MPa" else [A_min / (25.4 ** 2), A_max / (25.4 ** 2)]
#
#
#
#     -------------------------------------
#
#     # Import Standard Libraries
#     import logging
#
#
#     # ===========================================================================
#     #   Utilities
#     # ===========================================================================
#
#     def convert_2_psi(fck, units):
#         if units not in ["MPa", "psi"]:
#             raise KeyError("Only psi or MPa")
#         return fck if units == "psi" else 145.038 * fck
#
#
#     def convert_2_in(h, units):
#         if units not in ["MPa", "psi"]:
#             raise KeyError("Only psi or MPa")
#         return h if units == "psi" else h / 25.4
#
#
#     def convert_2_in2(A, units):
#         if units not in ["MPa", "psi"]:
#             raise KeyError("Only psi or MPa")
#         return A if units == "psi" else A / (25.4 ** 2)
#
#
#     def convert_2_pcf(rho, units):
#         if units not in ["MPa", "psi"]:
#             raise KeyError("Only psi or MPa")
#         if units == "psi":
#             return rho
#         elif units == "MPa" and rho != 145:
#             return 6.366 * rho
#
#
#     def convert_2_MPa(fck, units):
#         if units not in ["MPa", "psi"]:
#             raise KeyError("Only psi or MPa")
#         return fck if units == "MPa" else fck / 145.038
#
#
#     def convert_2_mm(h, units):
#         if units not in ["MPa", "psi"]:
#             raise KeyError("Only psi or MPa")
#         return h if units == "MPa" else 25.4 * h
#
#
#     def convert_2_mm2(A, units):
#         if units not in ["MPa", "psi"]:
#             raise KeyError("Only psi or MPa")
#         return A if units == "MPa" else 25.4 ** 2 * A
#
#
#     def compare_steel_area(As, A_min, A_max):
#         if As > A_min:
#             logging.info(
#                 "    As = {:6.2f} > A_min = {:6.2f}. Md > Mcr".format(As, A_min))
#         else:
#             logging.info(
#                 "    As = {:6.2f} < A_min = {:6.2f}. Md < Mcr".format(As, A_min))
#         if As < A_max:
#             logging.info(
#                 "    As = {:6.2f} < A_max = {:6.2f}. OK".format(As, A_max))
#         else:
#             logging.info(
#                 "    As = {:6.2f} > A_max = {:6.2f}. Not OK".format(As, A_max))