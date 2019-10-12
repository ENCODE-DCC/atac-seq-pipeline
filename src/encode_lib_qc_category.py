#!/usr/bin/env python3
"""
ENCODE QC log/plot to HTML converter

Author: Jin Lee (leepc12@gmail.com)
"""

from collections import OrderedDict
from base64 import b64encode


def to_number(var):
    """Convert to number or return None
    """
    try:
        if '.' in var:
            raise ValueError
        return int(var)
    except ValueError:
        try:
            return float(var)
        except ValueError:
            return None


class QCLog(object):
    """Parse a QC text file and convert it into a Python dict.

       TSV (number of columns >= 1) can be converted without a parser function.
       First column will be key name and the rest of columns will be
       values.

       For other QC log types, specify a parser function.
    """

    def __init__(self, log_file, parser=None):
        """
        Args:
            log_file: QC log file or QC dict
        """
        self._log_file = log_file
        self._parser = parser
        self._dict = None
        self.__parse()

    def to_dict(self):
        return self._dict

    def __parse(self):
        if isinstance(self._log_file, dict):
            self._dict = self._log_file
            return

        if self._parser is None:
            d = OrderedDict()
            # can parse TSV
            with open(self._log_file, 'r') as fp:
                lines = fp.read().strip('\n')

            def to_number_or_str(var):
                """Convert to number or return str
                """
                try:
                    if '.' in var:
                        raise ValueError
                    return int(var)
                except ValueError:
                    try:
                        return float(var)
                    except ValueError:
                        return var
                return None

            for i, line in enumerate(lines.split('\n')):
                arr = line.split('\t')
                if len(arr) == 1:
                    key = 'value' + str(i)
                    val = to_number(arr[0])
                elif len(arr) == 2:
                    key = arr[0]
                    val = to_number(arr[1])
                elif len(arr) > 2:
                    key = arr[0]
                    val = [to_number(v) for v in arr[1:]]
                else:
                    continue
                d[key] = val

            self._dict = d
        else:
            self._dict = self._parser(self._log_file)


class QCPlot(object):
    """Embed image as base64 string and return HTML string.
       QCPlot supports all images types supported by HTML <img>.
    """

    def __init__(self, plot_file, caption=None, size_pct=100):
        self._plot_file = plot_file
        self._caption = caption
        self._size_pct = size_pct
        self._encoded = None
        self._img_type = None
        self.__encode()

    def to_html(self):
        html = '''
        <figure style="display:inline-block">
          <img src="data:image/{img_type};base64,{encoded}" alt="{caption}" height="{size_pct}%"/>
          <figcaption style="text-align:center">{caption}</figcaption>
        </figure>
        '''
        return html.format(
            img_type=self._img_type,
            size_pct=self._size_pct,
            encoded=self._encoded,
            caption='' if self._caption is None else self._caption)

    def __encode(self):
        self._encoded = b64encode(
            open(self._plot_file, 'rb').read()).decode("utf-8")


class QCCategory(object):
    """QCCategory can have a child QCCategory and HTML will be resursively
    stacked. This is useful for having subcategories.
    """

    def __init__(self, cat_name, html_head='', html_foot='',
                 parser=None, map_key_desc=None, parent=None):
        """
        Args:
            cat_name: category name

            map_key_desc: map (keyname to description) for QCLog.
                          For example of samtools flagstat
                          'mapped_qc_failed' : 'Mapped(QC-failed)'

            parser: use it as default parser for all children QC logs
        """
        self._cat_name = cat_name
        self._html_head = html_head
        self._html_foot = html_foot
        self._parser = parser
        self._map_key_desc = map_key_desc
        self._qc_logs = OrderedDict()
        self._qc_plots = OrderedDict()
        self._child_categories = []
        if parent is not None:
            parent.add_category(self)

    def add_category(self, qc_category):
        self._child_categories.append(qc_category)

    def add_log(self, log_file, key=None):
        assert(key not in self._qc_logs)
        self._qc_logs[key] = QCLog(
            log_file,
            parser=self._parser,
        )

    def add_plot(self, plot_file, key=None, caption=None, size_pct=100):
        assert(key not in self._qc_plots)
        self._qc_plots[key] = QCPlot(
            plot_file,
            caption=caption if caption is not None else key,
            size_pct=size_pct
        )

    def to_dict(self):
        """Convert into dict. Plots will be ignored.
        """
        d = OrderedDict()
        for key, qc_log in self._qc_logs.items():
            d_ = qc_log.to_dict()
            if len(d_) > 0:
                d[key] = d_
        for cat in self._child_categories:
            d_ = cat.to_dict()
            if len(d_) > 0:
                d[cat._cat_name] = d_
        return d

    def to_html(self):
        """Print HTML only if there are contents to be shown
        """
        html = ''
        html += self.__qc_logs_to_html()
        html += self.__qc_plots_to_html()
        for cat in self._child_categories:
            html += cat.to_html()
        if html == '':
            return ''
        else:
            return self._html_head + html + self._html_foot

    def __qc_logs_to_html(self):
        """Print HTML only if there are contents to be shown
        Make an HTML table of qc_logs. For example,
                 rep1    rep2
        -------+-------+--------
        key1   | val1  | val1
        """
        if len(self._qc_logs) == 0:
            return ''

        html = '<table border="1" style="border-collapse:collapse">'

        # make HTML header row
        header = '<tr><th bgcolor="#EEEEEE">'
        arr = [' ']
        if len(self._qc_logs) == 1 and \
                list(self._qc_logs.keys())[0] is None:
            # skip header for single qc log
            arr += ['Description']
        else:
            arr += self._qc_logs.keys()
        header += '</th><th bgcolor="#EEEEEE">'.join(arr) + '</th></tr>\n'

        # declared as dict but will be used as set with empty values
        all_keys = OrderedDict()
        # contents
        for qc_log_k, qc_log_val in self._qc_logs.items():
            # print(qc_log_k, qc_log_val.to_dict())
            new_keys = OrderedDict.fromkeys(
                k for k in qc_log_val.to_dict() if k not in all_keys.keys())
            for new_key in new_keys:
                all_keys[new_key] = None

        content = ''
        for key in all_keys.keys():
            if self._map_key_desc is None:
                long_key_name = key
            else:
                long_key_name = self._map_key_desc[key]
            content += '<tr><th bgcolor="#EEEEEE" style="text-align:left">'\
                '{}</th><td>'.format(long_key_name)

            qc_log_content = []
            for qc_log_key, qc_log in self._qc_logs.items():
                d = qc_log.to_dict()
                if key not in d:
                    val = 'N/A'
                else:
                    val = str(d[key])
                qc_log_content.append(val)

            content += '</td><td>'.join(qc_log_content)
            content += '</td></tr>\n'

        html += header
        html += content
        html += '</table><br>\n'
        return html

    def __qc_plots_to_html(self):
        html = ''
        for k, qc_plot in self._qc_plots.items():
            html += qc_plot.to_html()
        return html
