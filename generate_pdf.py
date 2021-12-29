from fpdf import FPDF, HTMLMixin
from graphics_generation import graphics_generation
import spd_analyze
import os
import warnings
import copy
warnings.filterwarnings("ignore")

class CustomPDF(FPDF):

    def header(self):
        # Устанавливаем лого
        #self.image('pic_labels//koziy.png', 10, 2, 25)
        self.image('pic_labels//logo.png', 180, 6, 25)
        self.set_font('Arial', 'B', 15)

        self.cell(80)
        self.cell(0, 0, 'SPD Datasheet', ln=1)
        # Разрыв линии
        self.ln(20)

    def footer(self):
        self.set_y(-10)

        self.set_font('Arial', 'I', 8)

        # Добавляем номер страницы
        page = 'Page ' + str(self.page_no()) + '/{nb}'
        self.cell(0, 10, page, 0, 0, 'C')


class PdfDocument():
    def __init__(self, params):
        self.params = params
        self.SPD =  spd_analyze.SPD(self.params["dir_name"])
        self.SPD = self.SPD.get_spd()
        self.params["str_spd"] = ""

    def generate_pdf(self):
        self.criterion_list = ['SNR', 'QE']
        for criterion in self.criterion_list:
            params_specific_criterion = copy.deepcopy(self.params)
            params_specific_criterion["criterion"] = criterion
            self.generate_pdf_for_specific_criterion(params_specific_criterion)


    def generate_pdf_for_specific_criterion(self, params):
        self.pdf = CustomPDF()
        self.pdf.alias_nb_pages()

        if params["type"] == "gated":
            params["SPD"] = self.SPD
            grph = graphics_generation(params)
            grph.generate_pdf()
            params_spd = grph.params
            spd_pdf_class = GenerateGatedPdf(self.pdf, params_spd)
            spd_pdf_class.create_pdf()

        if params["type"] == "freerun":
            params["SPD"] = self.SPD
            grph = graphics_generation(params)
            grph.generate_pdf()
            params_spd = grph.params
            spd_pdf_class = GenerateFreerunPdf(self.pdf, params_spd)
            spd_pdf_class.create_pdf()

        if params["type"] == "butterfly":
            butterfly_params_dict = {}
            for key, spd in self.SPD.spd_dict.items():
                params_specific_spd = copy.deepcopy(params)
                type_spd = spd.type
                params_specific_spd["type"] = type_spd
                params_specific_spd["dir_name"] = os.path.join(params["dir_name"], self.SPD.path_dict[key])
                str_spd = "_" + str(key)
                params_specific_spd["str_spd"] = str_spd
                params_specific_spd["SPD"] = spd

                grph = graphics_generation(params_specific_spd)
                grph.generate_pdf()
                params_spd = grph.params
                butterfly_params_dict[key] = params_spd

                self.pdf = CustomPDF()
                self.pdf.alias_nb_pages()

                if type_spd == "gated":
                    spd_pdf_class = GenerateGatedPdf(self.pdf, params_specific_spd)
                if type_spd == "freerun":
                    spd_pdf_class = GenerateFreerunPdf(self.pdf, params_specific_spd)

                spd_pdf_class.create_pdf()

            self.pdf = CustomPDF()
            self.pdf.alias_nb_pages()

            spd_pdf_class = GenerateButterflyPdf(self.pdf, params, butterfly_params_dict, self.SPD)
            spd_pdf_class.create_pdf()


class GenerateButterflyPdf():
    def __init__(self, pdf, params, butterfly_params_dict,  SPD):
        self.pdf = pdf
        self.params = params
        self.butterfly_params_dict = butterfly_params_dict
        self.SPD = SPD

    def create_pdf(self):
        if self.params['criterion'] == 'QE':
            for key, spd in self.SPD.spd_dict.items():
                params_specific_spd = self.butterfly_params_dict[key]

                for el in params_specific_spd['temp']:
                    self.pdf.add_page()
                    self.pdf.set_font('Times', '', 14)

                    if params_specific_spd["type"] == "gated":
                        spd_pdf_class = GenerateGatedPdf(self.pdf, params_specific_spd)
                        filename = spd_pdf_class.create_list_criterion_qe_gated(el, str_spd=params_specific_spd["str_spd"])
                    elif params_specific_spd["type"] == "freerun":
                        spd_pdf_class = GenerateFreerunPdf(self.pdf, params_specific_spd)
                        filename = spd_pdf_class.create_list_criterion_qe_freerun(el, str_spd=params_specific_spd["str_spd"])

            path = os.path.join(self.params['dir_name'], params_specific_spd['fold'], self.params['criterion'])
            str_spd = ""
            if params_specific_spd['isvg'] == True:
                filename = os.path.join(path, 'report_QE_' + self.params['name'] + '_VG' + str_spd + '.pdf')
            else:
                filename = os.path.join(path, 'report_QE_' + self.params['name'] + '_CD' + str_spd + '.pdf')

            self.pdf.output(filename)

        if self.params['criterion'] == 'SNR':
            for key, spd in self.SPD.spd_dict.items():
                params_specific_spd = self.butterfly_params_dict[key]

                for el in params_specific_spd['temp']:
                    self.pdf.add_page()
                    self.pdf.set_font('Times', '', 14)

                    if params_specific_spd["type"] == "gated":
                        spd_pdf_class = GenerateGatedPdf(self.pdf, params_specific_spd)
                        filename = spd_pdf_class.create_list_criterion_snr_gated(el, str_spd=params_specific_spd[
                            "str_spd"])
                    elif params_specific_spd["type"] == "freerun":
                        spd_pdf_class = GenerateFreerunPdf(self.pdf, params_specific_spd)
                        filename = spd_pdf_class.create_list_criterion_snr_freerun(el, str_spd=params_specific_spd[
                            "str_spd"])

            path = os.path.join(self.params['dir_name'], params_specific_spd['fold'], self.params['criterion'])
            str_spd = ""
            if params_specific_spd['isvg'] == True:
                filename = os.path.join(path, 'report_QE_' + self.params['name'] + '_VG' + str_spd + '.pdf')
            else:
                filename = os.path.join(path, 'report_QE_' + self.params['name'] + '_CD' + str_spd + '.pdf')

            self.pdf.output(filename)



class GenerateGatedPdf():
    def __init__(self, pdf, params):
        self.pdf = pdf
        self.params = params

    def create_pdf(self):
        if self.params['criterion'] == 'QE':
            for el in self.params['temp']:
                self.pdf.add_page()
                self.pdf.set_font('Times', '', 14)

                filename = self.create_list_criterion_qe_gated(el, str_spd=self.params["str_spd"])

            self.pdf.output(filename)

        if self.params['criterion'] == 'SNR':
            for el in self.params['temp']:
                self.pdf.add_page()
                self.pdf.set_font('Times', '', 14)

                filename = self.create_list_criterion_snr_gated(el, str_spd=self.params["str_spd"])

            self.pdf.output(filename)

    def create_list_criterion_qe_gated(self, el, str_spd=""):
        w = 90
        self.pdf.cell(50, 0, txt="SPD type: " + self.params['type'], ln=1, align="C")
        self.pdf.cell(200, 0, txt="SPD Name: " + self.params['name'], ln=1, align="C")
        self.pdf.cell(350, 0, txt="Date: " + self.params['date'], ln=1, align="C")

        self.pdf.ln(10)
        self.pdf.cell(200, 0, txt="Optimal point", ln=1, align="C")
        self.pdf.cell(350, 0, txt="T = " + str(el) + ' C', ln=1, align="C")

        self.pdf.ln(10)
        simple_table(
            [self.params['optimal_params']['head'],
             self.params['optimal_params'][el][0],
             self.params['optimal_params'][el][1],
             self.params['optimal_params'][el][2]],
             self.pdf, criterion=self.params['criterion'])

        self.pdf.cell(30, 5, txt='', ln=1, align="C")
        # pdf.line(10, 68, 200, 68)
        simple_table([self.params['settings']['head'],
                      self.params['settings'][el][0],
                      self.params['settings'][el][1],
                      self.params['settings'][el][2]],
                     self.pdf,
                     criterion=self.params['criterion'])

        self.pdf.cell(200, 20, txt="Graphic characteristics representation", ln=1, align="C")

        if self.params['isvg'] == True:
            path = os.path.join(self.params['dir_name'], self.params['fold'], self.params['criterion'])
            self.pdf.image(path + '//HM_SNR_T' + str(el) + '_VG' + str_spd + '.png', x=10, y=135, w=w)
            # pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=120, w=90)
            self.pdf.image(path + '//HM_BIF_T' + str(el) + '_VG' + str_spd + '.png', x=110, y=135, w=w)
            self.pdf.image(path + '//DCR_T' + str(el) + '_VG' + str_spd + '.png', x=10, y=210, w=w)
            self.pdf.image(path + '//AP_T' + str(el) + '_VG' + str_spd + '.png', x=110, y=210, w=w)
        else:
            path = os.path.join(self.params['dir_name'], self.params['fold'], self.params['criterion'])
            self.pdf.image(path + '//HM_SNR_T' + str(el) + '_CD' + str_spd + '.png', x=10, y=135, w=w)
            # pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=120, w=90)
            self.pdf.image(path + '//HM_BIF_T' + str(el) + '_CD' + str_spd + '.png', x=110, y=135, w=w)
            self.pdf.image(path + '//DCR_T' + str(el) + '_CD' + str_spd + '.png', x=10, y=210, w=w)
            self.pdf.image(path + '//AP_T' + str(el) + '_CD' + str_spd + '.png', x=110, y=210, w=w)

        if self.params['isvg'] == True:
            filename = os.path.join(path, 'report_QE_' + self.params['name'] + '_VG' + str_spd + '.pdf')
        else:
            filename = os.path.join(path, 'report_QE_' + self.params['name'] + '_CD' + str_spd + '.pdf')

        return filename

    def create_list_criterion_snr_gated(self, el, str_spd=""):
        self.pdf.cell(50, 0, txt="SPD type: " + self.params['type'], ln=1, align="C")
        self.pdf.cell(200, 0, txt="SPD Name: " + self.params['name'], ln=1, align="C")
        self.pdf.cell(350, 0, txt="Date: " + self.params['date'], ln=1, align="C")

        self.pdf.ln(10)
        self.pdf.cell(200, 0, txt="Optimal point", ln=1, align="C")
        self.pdf.cell(350, 0, txt="T = " + str(el) + ' C', ln=1, align="C")

        self.pdf.ln(10)
        simple_table([self.params['optimal_params']['head'], self.params['optimal_params'][el]], self.pdf)
        self.pdf.cell(30, 5, txt='', ln=1, align="C")
        # pdf.line(10, 68, 200, 68)
        simple_table([self.params['settings']['head'], self.params['settings'][el]], self.pdf)

        self.pdf.cell(200, 20, txt="Graphic characteristics representation", ln=1, align="C")

        if self.params['isvg'] == True:
            path = os.path.join(self.params['dir_name'], self.params['fold'], self.params['criterion'])
            self.pdf.image(path + '//HM_SNR_T' + str(el) + '_VG' + str_spd + '.png', x=10, y=100, w=90)
            # pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=100, w=90)
            self.pdf.image(path + '//HM_BIF_T' + str(el) + '_VG' + str_spd + '.png', x=100, y=100, w=90)
            self.pdf.image(path + '//DCR_T' + str(el) + '_VG' + str_spd + '.png', x=10, y=180, w=90)
            self.pdf.image(path + '//AP_T' + str(el) + '_VG' + str_spd + '.png', x=100, y=180, w=90)
        else:
            path = os.path.join(self.params['dir_name'], self.params['fold'], self.params['criterion'])
            self.pdf.image(path + '//HM_SNR_T' + str(el) + '_CD' + str_spd + '.png', x=10, y=100, w=90)
            # pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=100, w=90)
            self.pdf.image(path + '//HM_BIF_T' + str(el) + '_CD' + str_spd + '.png', x=100, y=100, w=90)
            self.pdf.image(path + '//DCR_T' + str(el) + '_CD' + str_spd + '.png', x=10, y=180, w=90)
            self.pdf.image(path + '//AP_T' + str(el) + '_CD' + str_spd + '.png', x=100, y=180, w=90)

        if self.params['isvg'] == True:
            filename = os.path.join(path, 'report_SNR_' + self.params['name'] + '_VG' + str_spd + '.pdf')
        else:
            filename = os.path.join(path, 'report_SNR_' + self.params['name'] + '_CD' + str_spd + '.pdf')

        return filename

class GenerateFreerunPdf():
    def __init__(self, pdf, params):
        self.pdf = pdf
        self.params = params

    def create_pdf(self):
        if self.params['criterion'] == 'QE':
            for el in self.params['temp']:
                self.pdf.add_page()
                self.pdf.set_font('Times', '', 14)

                filename = self.create_list_criterion_qe_freerun(el, str_spd=self.params["str_spd"])

            self.pdf.output(filename)

        if self.params['criterion'] == 'SNR':
            for el in self.params['temp']:
                self.pdf.add_page()
                self.pdf.set_font('Times', '', 14)

                filename = self.create_list_criterion_snr_freerun(el, str_spd=self.params["str_spd"])

            self.pdf.output(filename)


    def create_list_criterion_qe_freerun(self, el, str_spd=""):
        self.pdf.cell(50, 0, txt="SPD type: " + self.params['type'], ln=1, align="C")
        self.pdf.cell(200, 0, txt="SPD Name: " + self.params['name'], ln=1, align="C")
        self.pdf.cell(350, 0, txt="Date: " + self.params['date'], ln=1, align="C")

        self.pdf.ln(10)
        self.pdf.cell(200, 0, txt="Optimal point", ln=1, align="C")
        self.pdf.cell(350, 0, txt="T = " + str(el) + ' C', ln=1, align="C")

        self.pdf.ln(10)
        simple_table([self.params['optimal_params']['head'], self.params['optimal_params'][el][0], self.params['optimal_params'][el][1]],
                     self.pdf, criterion=self.params['criterion'])
        self.pdf.cell(30, 5, txt='', ln=1, align="C")
        # pdf.line(10, 68, 200, 68)
        simple_table([self.params['settings']['head'], self.params['settings'][el][0], self.params['settings'][el][1]], self.pdf,
                     criterion=self.params['criterion'])

        self.pdf.cell(200, 20, txt="Graphic characteristics representation", ln=1, align="C")
        path = os.path.join(self.params['dir_name'], self.params['fold'], self.params['criterion'])
        self.pdf.image(path + '//SNR_T' + str(el) + str_spd + '.png', x=10, y=120, w=90)
        self.pdf.image(path + '//TR_T' + str(el) + str_spd + '.png', x=100, y=120, w=90)
        self.pdf.image(path + '//DCR_T' + str(el) + str_spd + '.png', x=10, y=200, w=95)
        self.pdf.image(path + '//AP_T' + str(el) + str_spd + '.png', x=110, y=200, w=85)
        filename = os.path.join(path, 'report_QE_' + self.params['name'] + str_spd + '.pdf')

        return filename


    def create_list_criterion_snr_freerun(self, el, str_spd=""):
        self.pdf.cell(50, 0, txt="SPD type: " + self.params['type'], ln=1, align="C")
        self.pdf.cell(200, 0, txt="SPD Name: " + self.params['name'], ln=1, align="C")
        self.pdf.cell(350, 0, txt="Date: " + self.params['date'], ln=1, align="C")

        self.pdf.ln(10)
        self.pdf.cell(200, 0, txt="Optimal point", ln=1, align="C")
        self.pdf.cell(350, 0, txt="T = " + str(el) + ' C', ln=1, align="C")

        self.pdf.ln(10)
        simple_table([self.params['optimal_params']['head'], self.params['optimal_params'][el]], self.pdf)
        self.pdf.cell(30, 5, txt='', ln=1, align="C")
        # pdf.line(10, 68, 200, 68)
        simple_table([self.params['settings']['head'], self.params['settings'][el]], self.pdf)

        self.pdf.cell(200, 20, txt="Graphic characteristics representation", ln=1, align="C")
        path = os.path.join(self.params['dir_name'], self.params['fold'], self.params['criterion'])
        self.pdf.image(path + '//SNR_T' + str(el) + str_spd + '.png', x=10, y=100, w=90)
        self.pdf.image(path + '//TR_T' + str(el) + str_spd + '.png', x=100, y=100, w=90)
        self.pdf.image(path + '//DCR_T' + str(el) + str_spd + '.png', x=10, y=180, w=95)
        self.pdf.image(path + '//AP_T' + str(el) + str_spd + '.png', x=110, y=180, w=85)
        filename = os.path.join(path, 'report_SNR_' + self.params['name'] + str_spd + '.pdf')

        return filename

def simple_table(data, pdf, spacing=1.5, criterion = 'SNR'):
    if criterion == 'QE':
        col_width = pdf.w / (len(data[0][0]) + 0.5)
        row_height = pdf.font_size
        for row in data:
            for item in row:
                pdf.cell(col_width, row_height * spacing,
                         txt=item, border=1)
            pdf.ln(row_height * spacing)

    if criterion == 'SNR':
        col_width = pdf.w / (len(data[0]) + 0.5)
        row_height = pdf.font_size
        for row in data:
            for item in row:
                pdf.cell(col_width, row_height * spacing,
                         txt=item, border=1)
            pdf.ln(row_height * spacing)

