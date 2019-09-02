#include "httpd.h"

#include "apr_strings.h"
#include "http_config.h"
#include "http_core.h"
#include "http_protocol.h"
#include "http_request.h"

static const char brunsli_mime_type[] = "image/x-j";

static int accept_brunsli(char* const accept) {
  static const char sep[] = ",";
  static const char qsep[] = ";";
  char* token;
  char* last;
  if (accept == NULL) return 0;
  for (token = apr_strtok(accept, sep, &last); token != NULL;
       token = apr_strtok(NULL, sep, &last)) {
    char* qlast;
    char* const mime_type = apr_strtok(token, qsep, &qlast);
    apr_collapse_spaces(mime_type, mime_type);
    if (strcmp(mime_type, brunsli_mime_type) == 0) {
      return 1;
    }
  }
  return 0;
}

static int brunsli_handler(request_rec* const r) {
  int rc;
  apr_finfo_t finfo;
  char* filename;
  char* const accept =
      apr_pstrdup(r->pool, apr_table_get(r->headers_in, "Accept"));
  if (!accept_brunsli(accept)) {
    return DECLINED;
  }

  rc = ap_core_translate(r);
  if (rc != APR_SUCCESS || r->filename == NULL) {
    return DECLINED;
  }

  filename = apr_psprintf(r->pool, "%s.j", r->filename);

  rc = apr_stat(&finfo, filename, APR_FINFO_MIN, r->pool);
  if (!(rc == APR_SUCCESS && finfo.filetype != APR_NOFILE &&
        (finfo.filetype & APR_DIR) == 0)) {
    return DECLINED;
  }

  r->filename = filename;
  ap_set_content_type(r, brunsli_mime_type);

  return OK;
}

static void register_hooks(apr_pool_t* const pool) {
  static const char* const predecessors[] = {"http_core.c", NULL};
  ap_hook_translate_name(brunsli_handler, predecessors, NULL, APR_HOOK_MIDDLE);
}

AP_DECLARE_MODULE(brunsli) = {
    STANDARD20_MODULE_STUFF, NULL, NULL, NULL, NULL, NULL, register_hooks};
