from fpdf import FPDF, HTMLMixin
import os

class CustomPDF(FPDF):

    def header(self):
        # Устанавливаем лого
        self.image('pic_labels//koziy.png', 10, 2, 25)
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

    #First page
    pdf.add_page()
    pdf.set_font('Times', '', 14)

    pdf.cell(200, 0, txt="SPD Name: " + params['name'] , ln=1, align="C")
    pdf.cell(350, 0, txt="Date: " + params['date'], ln=1, align="C")

    pdf.cell(200, 20, txt="Optimal point", ln=1, align="C")

    simple_table(params['optimal_params'], pdf)
    pdf.cell(30,5, txt='', ln=1, align="C")
    #pdf.line(10, 68, 200, 68)
    simple_table(params['settings'], pdf)

    pdf.cell(200, 20, txt="DCR dependence on gate", ln=1, align="C")

    path = os.path.join(params['dir_name'], params['fold'])
    pdf.image(path + '//DCR_PDE10.png', x=10, y=100, w=90)
    pdf.image(path + '//DCR_PDE15.png', x=100, y=100, w=90)
    pdf.image(path + '//DCR_PDE20.png', x=10, y=180, w=90)
    pdf.image(path + '//DCR_PDE25.png', x=100, y=180, w=90)

    #Second page
    pdf.add_page()
    pdf.cell(200, 0, txt="SPD Name: " + params['name'], ln=1, align="C")
    pdf.cell(350, 0, txt="Date: " + params['date'], ln=1, align="C")

    pdf.cell(200, 20, txt="AP dependence on gate", ln=1, align="C")

    pdf.image(path + '//AP_PDE10.png', x=10, y=50, w=90)
    pdf.image(path + '//AP_PDE15.png', x=100, y=50, w=90)
    pdf.image(path + '//AP_PDE20.png', x=10, y=125, w=90)
    pdf.image(path + '//AP_PDE25.png', x=100, y=125, w=90)

    pdf.ln(155)
    pdf.cell(100, 0, txt="3D map", ln=1, align="C")
    pdf.cell(280, 0, txt="Optimal pont TR", ln=1, align="C")

    pdf.image(path + '//HM_SNR_T-55.png', x=10, y=210, w=90)
    pdf.image(path + '//TR.png', x=100, y=210, w=90)


    pdf.output(os.path.join(path,'report_' + params['name'] + '.pdf' ))


def simple_table(data, pdf, spacing=1.5):
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