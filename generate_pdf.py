from fpdf import FPDF, HTMLMixin
import os

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


def create_pdf(params):
    pdf = CustomPDF()
    # Создаем особое значение {nb}
    pdf.alias_nb_pages()

    if params['criterion'] == 'QE':
        for el in params['temp']:
            pdf.add_page()
            pdf.set_font('Times', '', 14)

            if params['type'] == 'gated':
                pdf.cell(50, 0, txt="SPD type: " + params['type'], ln=1, align="C")
                pdf.cell(200, 0, txt="SPD Name: " + params['name'], ln=1, align="C")
                pdf.cell(350, 0, txt="Date: " + params['date'], ln=1, align="C")

                pdf.ln(10)
                pdf.cell(200, 0, txt="Optimal point", ln=1, align="C")
                pdf.cell(350, 0, txt="T = " + str(el) + ' C', ln=1, align="C")

                pdf.ln(10)
                simple_table([params['optimal_params']['head'], params['optimal_params'][el][0], params['optimal_params'][el][1]], pdf, criterion=params['criterion'])
                pdf.cell(30, 5, txt='', ln=1, align="C")
                # pdf.line(10, 68, 200, 68)
                simple_table([params['settings']['head'], params['settings'][el][0], params['settings'][el][1]], pdf, criterion=params['criterion'])

                pdf.cell(200, 20, txt="Graphic characteristics representation", ln=1, align="C")

                if params['isvg'] == True:
                    path = os.path.join(params['dir_name'], params['fold'], params['criterion'])
                    pdf.image(path + '//HM_SNR_T' + str(el) + '_VG' + '.png', x=10, y=120, w=90)
                    #pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=120, w=90)
                    pdf.image(path + '//HM_BIF_T' + str(el) + '_VG' +  '.png', x=100, y=120, w=90)
                    pdf.image(path + '//DCR_T' + str(el)  + '_VG' +  '.png', x=10, y=200, w=90)
                    pdf.image(path + '//AP_T' + str(el)  + '_VG' +  '.png', x=100, y=200, w=90)
                else:
                    path = os.path.join(params['dir_name'], params['fold'], params['criterion'])
                    pdf.image(path + '//HM_SNR_T' + str(el)  + '_CD' +  '.png', x=10, y=120, w=90)
                    # pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=120, w=90)
                    pdf.image(path + '//HM_BIF_T' + str(el) + '_CD' + '.png', x=100, y=120, w=90)
                    pdf.image(path + '//DCR_T' + str(el) + '_CD' + '.png', x=10, y=200, w=90)
                    pdf.image(path + '//AP_T' + str(el) + '_CD' + '.png', x=100, y=200, w=90)

                if params['isvg'] == True:
                    filename = os.path.join(path, 'report_QE_' + params['name'] + '_VG' + '.pdf')
                else:
                    filename = os.path.join(path, 'report_QE_' + params['name'] + '_CD' + '.pdf')

            elif params['type'] == 'freerun':
                pdf.cell(50, 0, txt="SPD type: " + params['type'], ln=1, align="C")
                pdf.cell(200, 0, txt="SPD Name: " + params['name'], ln=1, align="C")
                pdf.cell(350, 0, txt="Date: " + params['date'], ln=1, align="C")

                pdf.ln(10)
                pdf.cell(200, 0, txt="Optimal point", ln=1, align="C")
                pdf.cell(350, 0, txt="T = " + str(el) + ' C', ln=1, align="C")

                pdf.ln(10)
                simple_table([params['optimal_params']['head'], params['optimal_params'][el][0], params['optimal_params'][el][1]], pdf, criterion=params['criterion'])
                pdf.cell(30, 5, txt='', ln=1, align="C")
                # pdf.line(10, 68, 200, 68)
                simple_table([params['settings']['head'], params['settings'][el][0], params['settings'][el][1]], pdf, criterion=params['criterion'])

                pdf.cell(200, 20, txt="Graphic characteristics representation", ln=1, align="C")
                path = os.path.join(params['dir_name'], params['fold'], params['criterion'])
                pdf.image(path + '//SNR_T' + str(el) + '.png', x=10, y=120, w=90)
                pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=120, w=90)
                pdf.image(path + '//DCR_T' + str(el) + '.png', x=10, y=200, w=95)
                pdf.image(path + '//AP_T' + str(el) + '.png', x=110, y=200, w=85)
                filename = os.path.join(path, 'report_QE_' + params['name'] + '.pdf')

        pdf.output(filename)

    if params['criterion'] == 'SNR':
        for el in params['temp']:
            pdf.add_page()
            pdf.set_font('Times', '', 14)

            if params['type'] == 'gated':
                pdf.cell(50, 0, txt="SPD type: " + params['type'], ln=1, align="C")
                pdf.cell(200, 0, txt="SPD Name: " + params['name'] , ln=1, align="C")
                pdf.cell(350, 0, txt="Date: " + params['date'], ln=1, align="C")

                pdf.ln(10)
                pdf.cell(200, 0, txt="Optimal point", ln=1, align="C")
                pdf.cell(350, 0, txt="T = " + str(el) + ' C', ln=1, align="C")

                pdf.ln(10)
                simple_table([params['optimal_params']['head'],params['optimal_params'][el]], pdf)
                pdf.cell(30,5, txt='', ln=1, align="C")
                #pdf.line(10, 68, 200, 68)
                simple_table([params['settings']['head'], params['settings'][el]], pdf)

                pdf.cell(200, 20, txt="Graphic characteristics representation", ln=1, align="C")

                if params['isvg'] == True:
                    path = os.path.join(params['dir_name'], params['fold'], params['criterion'])
                    pdf.image(path + '//HM_SNR_T' + str(el) + '_VG' + '.png', x=10, y=100, w=90)
                    #pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=100, w=90)
                    pdf.image(path + '//HM_BIF_T' + str(el) + '_VG' + '.png', x=100, y=100, w=90)
                    pdf.image(path + '//DCR_T' + str(el) + '_VG' + '.png', x=10, y=180, w=90)
                    pdf.image(path + '//AP_T' + str(el) + '_VG' + '.png', x=100, y=180, w=90)
                else:
                    path = os.path.join(params['dir_name'], params['fold'], params['criterion'])
                    pdf.image(path + '//HM_SNR_T' + str(el) + '_CD' + '.png', x=10, y=100, w=90)
                    # pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=100, w=90)
                    pdf.image(path + '//HM_BIF_T' + str(el) + '_CD' + '.png', x=100, y=100, w=90)
                    pdf.image(path + '//DCR_T' + str(el) + '_CD' + '.png', x=10, y=180, w=90)
                    pdf.image(path + '//AP_T' + str(el) + '_CD' + '.png', x=100, y=180, w=90)

                if params['isvg'] == True:
                    filename = os.path.join(path, 'report_SNR_' + params['name'] + '_VG' + '.pdf')
                else:
                    filename = os.path.join(path, 'report_SNR_' + params['name'] + '_CD' + '.pdf')

            elif params['type'] == 'freerun':
                pdf.cell(50, 0, txt="SPD type: " + params['type'], ln=1, align="C")
                pdf.cell(200, 0, txt="SPD Name: " + params['name'], ln=1, align="C")
                pdf.cell(350, 0, txt="Date: " + params['date'], ln=1, align="C")

                pdf.ln(10)
                pdf.cell(200, 0, txt="Optimal point", ln=1, align="C")
                pdf.cell(350, 0, txt="T = " + str(el) + ' C', ln=1, align="C")

                pdf.ln(10)
                simple_table([params['optimal_params']['head'], params['optimal_params'][el]], pdf)
                pdf.cell(30, 5, txt='', ln=1, align="C")
                # pdf.line(10, 68, 200, 68)
                simple_table([params['settings']['head'], params['settings'][el]], pdf)

                pdf.cell(200, 20, txt="Graphic characteristics representation", ln=1, align="C")
                path = os.path.join(params['dir_name'], params['fold'], params['criterion'])
                pdf.image(path + '//SNR_T' + str(el) + '.png', x=10, y=100, w=90)
                pdf.image(path + '//TR_T' + str(el) + '.png', x=100, y=100, w=90)
                pdf.image(path + '//DCR_T' + str(el) + '.png', x=10, y=180, w=95)
                pdf.image(path + '//AP_T' + str(el) + '.png', x=110, y=180, w=85)
                filename = os.path.join(path, 'report_SNR_' + params['name'] + '.pdf')

        pdf.output(filename)


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




'''
params = {}
params['name'] = 'Uliana'
params['date'] = '05.06.2021'
params['optimal_params'] = [['Criteria', 'PDE, %', 'DCR, Hz', 'SNR', 'AP, %', 'DT, mus', 'TR, ps'],
            ['SNR', '20 %', '200', '0.18', '3', '5', '200']]
params['settings'] = [['Criteria','Vg', 'Vb', 'T'], ['SNR', '3ff', '64.55', '-55']]
params['temp_cart'] = 'temp_cart.png'
params['timeres'] = 'timeres.png'
params['dcr(pde)'] = 'dcr(pde).png'
params['ap(pde)'] = 'ap(pde).png'

create_pdf('my_pdf_1.pdf', params)
'''